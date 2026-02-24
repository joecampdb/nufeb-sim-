"""
Mercury Reef — NUFEB VTU to Blender Importer
=============================================
Imports NUFEB-2 VTK particle data into Blender for photorealistic
Cycles rendering of biofilm simulations.

Usage:
  1. Open Blender 4.0+
  2. Switch to the Scripting workspace
  3. Open this file (or paste it)
  4. Edit the CONFIGURATION section below (set your VTU path)
  5. Click "Run Script" (or Alt+P)

The script will:
  - Parse the VTU file for particle positions, diameters, and types
  - Create one Blender object per species type (HET, AOB, NOB, etc.)
  - Set up Geometry Nodes on each for sphere instancing (scaled by diameter)
  - Create Principled BSDF materials with subsurface scattering
  - Optionally set up camera, lighting, and substrate plane

Author: Claude Code / Anthropic
Project: Hg(SCN)2 Bioremediation by Anammox Biofilm
"""

import bpy
import bmesh
import os
import xml.etree.ElementTree as ET
import numpy as np
from mathutils import Vector

# ╔══════════════════════════════════════════════════════════════╗
# ║                     CONFIGURATION                           ║
# ╚══════════════════════════════════════════════════════════════╝

# Path to a single VTU file (use raw string on Windows)
VTU_PATH = r"C:\Users\ltjjp\OneDrive\Desktop\mercury_reef_vtk_combined\dump3760.vtu"

# Scale factor: NUFEB uses SI meters; multiply by 1e6 so
# 1 Blender unit = 1 micrometer. The 400um domain becomes
# 400 Blender units -- a comfortable working size.
SCALE = 1e6

# Sphere resolution for instancing (faces per sphere).
# Lower = faster viewport, higher = smoother render.
#   1 = 20 faces  (draft/preview)
#   2 = 80 faces  (good quality -- recommended)
#   3 = 320 faces (high quality, slow with >500k particles)
SPHERE_SUBDIVISIONS = 2

# Dead cell decimation: only import every Nth dead cell.
# Set to 1 to import all dead cells (slow for large frames).
# Set to 10 to import 10% of dead cells (fast, still shows necrotic core).
# Set to 0 to skip dead cells entirely.
DEAD_DECIMATE = 5

# Whether to set up camera, lights, and substrate plane
SETUP_SCENE = True

# Whether to clear the existing scene before importing
CLEAR_SCENE = True


# ╔══════════════════════════════════════════════════════════════╗
# ║                   SPECIES DEFINITIONS                       ║
# ╚══════════════════════════════════════════════════════════════╝

# Each species: (name, base_color_rgb, shader_properties)
# Colors match the earth-tone palette from the NUFEB simulation.
# Shader properties map to Principled BSDF inputs in Blender 4.0.

SPECIES = {
    1: {
        "name": "HET",
        "color": (0.824, 0.706, 0.549),  # Tan
        "shader": {
            "Subsurface Weight": 0.1,
            "Subsurface Radius": (0.8, 0.5, 0.3),
            "Roughness": 0.4,
            "Coat Weight": 0.15,
            "Coat Roughness": 0.1,
        },
    },
    2: {
        "name": "AOB",
        "color": (0.333, 0.420, 0.184),  # Dark Olive Green
        "shader": {
            "Subsurface Weight": 0.15,
            "Subsurface Radius": (0.4, 0.7, 0.3),
            "Roughness": 0.35,
            "Coat Weight": 0.2,
            "Coat Roughness": 0.1,
        },
    },
    3: {
        "name": "NOB",
        "color": (0.627, 0.322, 0.176),  # Sienna
        "shader": {
            "Subsurface Weight": 0.1,
            "Subsurface Radius": (0.7, 0.4, 0.2),
            "Roughness": 0.4,
            "Coat Weight": 0.15,
            "Coat Roughness": 0.1,
        },
    },
    4: {
        "name": "ANA",
        "color": (1.000, 0.843, 0.000),  # Gold
        "shader": {
            "Subsurface Weight": 0.2,
            "Subsurface Radius": (1.0, 0.6, 0.1),
            "Roughness": 0.3,
            "Coat Weight": 0.25,
            "Coat Roughness": 0.05,
        },
    },
    5: {
        "name": "CECT",
        "color": (0.000, 0.502, 0.502),  # Teal
        "shader": {
            "Subsurface Weight": 0.15,
            "Subsurface Radius": (0.2, 0.6, 0.7),
            "Roughness": 0.35,
            "Coat Weight": 0.2,
            "Coat Roughness": 0.1,
        },
    },
    6: {
        "name": "EPS",
        "color": (0.961, 0.871, 0.702),  # Wheat
        "shader": {
            # Full subsurface + transmission: translucent wet gel
            "Subsurface Weight": 1.0,
            "Subsurface Radius": (1.5, 1.0, 0.5),
            "Subsurface Scale": 0.5,
            "Roughness": 0.15,
            "Transmission Weight": 0.4,
            "IOR": 1.33,
            "Alpha": 0.8,
        },
    },
    7: {
        "name": "DEAD",
        "color": (0.412, 0.412, 0.412),  # Dark Gray
        "shader": {
            # Matte, opaque, lifeless
            "Subsurface Weight": 0.0,
            "Roughness": 0.85,
            "Specular IOR Level": 0.3,
        },
    },
}


# ╔══════════════════════════════════════════════════════════════╗
# ║                    VTU FILE PARSER                          ║
# ╚══════════════════════════════════════════════════════════════╝

def parse_vtu(filepath):
    """
    Parse a LAMMPS/NUFEB ASCII VTU file.

    Returns:
        positions : np.ndarray (N, 3) -- particle xyz coordinates
        types     : np.ndarray (N,)   -- integer species type
        diameters : np.ndarray (N,)   -- particle diameters
    """
    print(f"Parsing VTU: {filepath}")
    print(f"  File size: {os.path.getsize(filepath) / 1e6:.1f} MB")

    tree = ET.parse(filepath)
    root = tree.getroot()

    piece = root.find(".//Piece")
    n_points = int(piece.get("NumberOfPoints"))
    print(f"  Particles: {n_points:,}")

    # Extract PointData arrays
    point_data = piece.find("PointData")
    data_arrays = {}
    for da in point_data.findall("DataArray"):
        name = da.get("Name")
        dtype = da.get("type")
        text = da.text.strip()

        if dtype in ("Float64", "Float32"):
            data_arrays[name] = np.fromstring(text, sep=" ", dtype=np.float64)
        elif dtype in ("Int32", "Int64"):
            data_arrays[name] = np.fromstring(text, sep=" ", dtype=np.int32)

    # Extract point coordinates
    points_elem = piece.find(".//Points/DataArray")
    n_components = int(points_elem.get("NumberOfComponents", 3))
    positions_flat = np.fromstring(points_elem.text.strip(), sep=" ", dtype=np.float64)
    positions = positions_flat.reshape(-1, n_components)

    types = data_arrays["type"]
    diameters = data_arrays["diameter"]

    print(f"  Type range: {types.min()} - {types.max()}")
    print(f"  Diameter range: {diameters.min():.2e} - {diameters.max():.2e} m")

    # Report per-type counts
    for t in sorted(np.unique(types)):
        count = np.sum(types == t)
        name = SPECIES.get(t, {}).get("name", f"Type{t}")
        print(f"    {name} (type {t}): {count:,}")

    return positions, types, diameters


# ╔══════════════════════════════════════════════════════════════╗
# ║                  MATERIAL CREATION                          ║
# ╚══════════════════════════════════════════════════════════════╝

def create_material(species_id):
    """Create a Principled BSDF material for a species type."""
    spec = SPECIES[species_id]
    mat_name = f"MAT_{spec['name']}"

    # Reuse existing material if it exists
    if mat_name in bpy.data.materials:
        return bpy.data.materials[mat_name]

    mat = bpy.data.materials.new(name=mat_name)
    mat.use_nodes = True

    # Get the Principled BSDF node
    bsdf = mat.node_tree.nodes.get("Principled BSDF")
    if bsdf is None:
        bsdf = mat.node_tree.nodes.new("ShaderNodeBsdfPrincipled")

    # Set base color
    r, g, b = spec["color"]
    bsdf.inputs["Base Color"].default_value = (r, g, b, 1.0)

    # Set shader properties
    for prop, value in spec["shader"].items():
        if prop in bsdf.inputs:
            inp = bsdf.inputs[prop]
            if isinstance(value, tuple):
                # Vector input (e.g., Subsurface Radius)
                inp.default_value = value
            else:
                inp.default_value = value

    # For EPS: enable alpha blending
    if species_id == 6:
        mat.blend_method = 'BLEND'  # Eevee transparency
        mat.use_backface_culling = False

    return mat


# ╔══════════════════════════════════════════════════════════════╗
# ║              BLENDER OBJECT CREATION                        ║
# ╚══════════════════════════════════════════════════════════════╝

def create_point_cloud(name, positions, radii):
    """
    Create a Blender mesh object with vertices at particle positions
    and a 'radius' float attribute for Geometry Nodes scaling.
    """
    mesh = bpy.data.meshes.new(f"Mesh_{name}")
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)

    # Build mesh from vertices only (no edges or faces)
    bm = bmesh.new()
    for pos in positions:
        bm.verts.new(pos)
    bm.to_mesh(mesh)
    bm.free()

    # Store radius as a vertex attribute for Geometry Nodes
    attr = mesh.attributes.new(name="radius", type='FLOAT', domain='POINT')
    attr.data.foreach_set("value", radii.tolist())

    return obj


def create_sphere_template():
    """Create a small ico sphere mesh to use as the instance template."""
    name = "_SphereTemplate"
    if name in bpy.data.objects:
        return bpy.data.objects[name]

    mesh = bpy.data.meshes.new(f"Mesh{name}")
    obj = bpy.data.objects.new(name, mesh)
    bpy.context.collection.objects.link(obj)

    bm = bmesh.new()
    bmesh.ops.create_icosphere(
        bm,
        subdivisions=SPHERE_SUBDIVISIONS,
        radius=1.0,  # Unit radius; Geometry Nodes scales to actual size
    )
    bm.to_mesh(mesh)
    bm.free()

    # Smooth shading
    for poly in mesh.polygons:
        poly.use_smooth = True

    # Hide template from viewport and render
    obj.hide_set(True)
    obj.hide_render = True

    return obj


# ╔══════════════════════════════════════════════════════════════╗
# ║              GEOMETRY NODES SETUP                           ║
# ╚══════════════════════════════════════════════════════════════╝

def setup_geometry_nodes(obj, sphere_obj, material):
    """
    Create a Geometry Nodes modifier that instances spheres at each
    vertex, scaled by the 'radius' attribute, with the given material.
    """
    # Create node group
    group_name = f"GN_{obj.name}"
    if group_name in bpy.data.node_groups:
        bpy.data.node_groups.remove(bpy.data.node_groups[group_name])

    group = bpy.data.node_groups.new(group_name, 'GeometryNodeTree')

    # Create interface sockets (Blender 4.0+)
    group.interface.new_socket(
        name="Geometry", in_out='INPUT', socket_type='NodeSocketGeometry'
    )
    group.interface.new_socket(
        name="Geometry", in_out='OUTPUT', socket_type='NodeSocketGeometry'
    )

    # Clear default nodes
    group.nodes.clear()

    # -- Create nodes --

    # Group Input
    n_input = group.nodes.new('NodeGroupInput')
    n_input.location = (-800, 0)

    # Group Output
    n_output = group.nodes.new('NodeGroupOutput')
    n_output.location = (600, 0)

    # Object Info (references the sphere template)
    n_obj_info = group.nodes.new('GeometryNodeObjectInfo')
    n_obj_info.location = (-500, -200)
    n_obj_info.inputs['Object'].default_value = sphere_obj
    n_obj_info.transform_space = 'RELATIVE'

    # Named Attribute: read "radius" from vertex data
    n_radius = group.nodes.new('GeometryNodeInputNamedAttribute')
    n_radius.location = (-500, -50)
    n_radius.data_type = 'FLOAT'
    n_radius.inputs['Name'].default_value = "radius"

    # Combine XYZ: convert scalar radius to (r, r, r) vector for scaling
    n_combine = group.nodes.new('ShaderNodeCombineXYZ')
    n_combine.location = (-300, -50)

    # Instance on Points
    n_instance = group.nodes.new('GeometryNodeInstanceOnPoints')
    n_instance.location = (-100, 0)

    # Set Material
    n_set_mat = group.nodes.new('GeometryNodeSetMaterial')
    n_set_mat.location = (200, 0)
    n_set_mat.inputs['Material'].default_value = material

    # -- Link nodes --
    links = group.links

    # Geometry flow: Input -> Instance on Points -> Set Material -> Output
    links.new(n_input.outputs[0], n_instance.inputs['Points'])
    links.new(n_instance.outputs['Instances'], n_set_mat.inputs['Geometry'])
    links.new(n_set_mat.outputs['Geometry'], n_output.inputs[0])

    # Sphere template: Object Info -> Instance on Points
    links.new(n_obj_info.outputs['Geometry'], n_instance.inputs['Instance'])

    # Radius -> Combine XYZ -> Scale
    links.new(n_radius.outputs['Attribute'], n_combine.inputs['X'])
    links.new(n_radius.outputs['Attribute'], n_combine.inputs['Y'])
    links.new(n_radius.outputs['Attribute'], n_combine.inputs['Z'])
    links.new(n_combine.outputs['Vector'], n_instance.inputs['Scale'])

    # Apply modifier to object
    mod = obj.modifiers.new("GeometryNodes", 'NODES')
    mod.node_group = group

    return group


# ╔══════════════════════════════════════════════════════════════╗
# ║                    SCENE SETUP                              ║
# ╚══════════════════════════════════════════════════════════════╝

def setup_scene(domain_center, domain_size):
    """Set up camera, lighting, and substrate plane."""

    cx, cy, cz = domain_center
    dx, dy, dz = domain_size

    # --- Substrate plane ---
    bpy.ops.mesh.primitive_plane_add(
        size=max(dx, dy) * 2,
        location=(cx, cy, 0),
    )
    plane = bpy.context.active_object
    plane.name = "Substrate"

    # Substrate material: dark matte surface
    mat = bpy.data.materials.new("MAT_Substrate")
    mat.use_nodes = True
    bsdf = mat.node_tree.nodes["Principled BSDF"]
    bsdf.inputs["Base Color"].default_value = (0.15, 0.13, 0.12, 1.0)
    bsdf.inputs["Roughness"].default_value = 0.7
    plane.data.materials.append(mat)

    # --- Camera ---
    cam_data = bpy.data.cameras.new("Camera")
    cam_data.lens = 50
    cam_data.dof.use_dof = True
    cam_data.dof.aperture_fstop = 4.0

    cam_obj = bpy.data.objects.new("Camera", cam_data)
    bpy.context.collection.objects.link(cam_obj)

    # Position camera at 45 degrees, looking at biofilm center
    biofilm_top = cz * 0.5  # Approximate biofilm center height
    cam_dist = max(dx, dy) * 2.0
    cam_obj.location = (
        cx + cam_dist * 0.6,
        cy - cam_dist * 0.8,
        biofilm_top + cam_dist * 0.5,
    )

    # Point camera at biofilm center
    direction = Vector((cx, cy, biofilm_top)) - cam_obj.location
    rot_quat = direction.to_track_quat('-Z', 'Y')
    cam_obj.rotation_euler = rot_quat.to_euler()

    # Set focus distance to biofilm center
    cam_data.dof.focus_distance = direction.length

    bpy.context.scene.camera = cam_obj

    # --- Key Light ---
    key_data = bpy.data.lights.new("KeyLight", 'AREA')
    key_data.energy = 500
    key_data.size = max(dx, dy)
    key_data.color = (1.0, 0.98, 0.95)

    key_obj = bpy.data.objects.new("KeyLight", key_data)
    bpy.context.collection.objects.link(key_obj)
    key_obj.location = (cx + dx, cy - dy * 0.5, dz * 1.5)
    key_obj.rotation_euler = (0.8, 0, 0.4)

    # --- Fill Light ---
    fill_data = bpy.data.lights.new("FillLight", 'AREA')
    fill_data.energy = 200
    fill_data.size = max(dx, dy) * 0.8
    fill_data.color = (0.85, 0.9, 1.0)

    fill_obj = bpy.data.objects.new("FillLight", fill_data)
    bpy.context.collection.objects.link(fill_obj)
    fill_obj.location = (cx - dx * 0.8, cy + dy, dz)
    fill_obj.rotation_euler = (1.0, 0, -0.6)

    # --- Rim Light ---
    rim_data = bpy.data.lights.new("RimLight", 'AREA')
    rim_data.energy = 300
    rim_data.size = max(dx, dy) * 0.5

    rim_obj = bpy.data.objects.new("RimLight", rim_data)
    bpy.context.collection.objects.link(rim_obj)
    rim_obj.location = (cx - dx * 0.3, cy + dy * 0.5, dz * 2)
    rim_obj.rotation_euler = (0.3, 0, -0.2)

    # --- Render settings ---
    scene = bpy.context.scene
    scene.render.engine = 'CYCLES'
    scene.cycles.samples = 256          # Bump to 1024+ for final render
    scene.cycles.use_denoising = True
    scene.render.resolution_x = 1920
    scene.render.resolution_y = 1080
    scene.render.film_transparent = True  # Transparent background

    # Light bounces for SSS and transmission (EPS needs these)
    scene.cycles.max_bounces = 12
    scene.cycles.diffuse_bounces = 4
    scene.cycles.glossy_bounces = 4
    scene.cycles.transmission_bounces = 8
    scene.cycles.transparent_max_bounces = 8

    # World: subtle gray environment (swap for HDRI for production)
    world = bpy.data.worlds.get("World")
    if world is None:
        world = bpy.data.worlds.new("World")
    scene.world = world
    world.use_nodes = True
    bg = world.node_tree.nodes.get("Background")
    if bg:
        bg.inputs["Color"].default_value = (0.05, 0.05, 0.06, 1.0)
        bg.inputs["Strength"].default_value = 0.3

    print("Scene setup complete (camera, 3-point lighting, Cycles render)")
    print("  TIP: For production renders, add an HDRI from polyhaven.com")
    print("       to the World shader for natural environment lighting.")


# ╔══════════════════════════════════════════════════════════════╗
# ║                      MAIN                                   ║
# ╚══════════════════════════════════════════════════════════════╝

def main():
    print("\n" + "=" * 60)
    print("  Mercury Reef -- NUFEB VTU to Blender Importer")
    print("=" * 60 + "\n")

    # --- Clear scene ---
    if CLEAR_SCENE:
        bpy.ops.object.select_all(action='SELECT')
        bpy.ops.object.delete(use_global=True)
        # Clean orphan data
        for mesh in bpy.data.meshes:
            bpy.data.meshes.remove(mesh)
        for mat in bpy.data.materials:
            bpy.data.materials.remove(mat)
        for ng in bpy.data.node_groups:
            bpy.data.node_groups.remove(ng)
        print("Scene cleared.\n")

    # --- Parse VTU ---
    positions, types, diameters = parse_vtu(VTU_PATH)

    # Scale to micrometers
    positions_scaled = positions * SCALE
    radii_scaled = (diameters * SCALE) / 2.0  # Diameter -> radius

    # Domain bounds (for scene setup)
    pos_min = positions_scaled.min(axis=0)
    pos_max = positions_scaled.max(axis=0)
    domain_center = (pos_min + pos_max) / 2.0
    domain_size = pos_max - pos_min

    print(f"\nDomain (um): {domain_size[0]:.0f} x {domain_size[1]:.0f} x {domain_size[2]:.0f}")
    print(f"Center (um): ({domain_center[0]:.0f}, {domain_center[1]:.0f}, {domain_center[2]:.0f})\n")

    # --- Create sphere template ---
    sphere_obj = create_sphere_template()
    print(f"Sphere template: {SPHERE_SUBDIVISIONS} subdivisions "
          f"({len(sphere_obj.data.polygons)} faces)\n")

    # --- Create per-species objects ---
    created_objects = []
    total_imported = 0

    for type_id in sorted(SPECIES.keys()):
        spec = SPECIES[type_id]
        mask = types == type_id
        count = mask.sum()

        if count == 0:
            print(f"  {spec['name']}: 0 particles -- skipping")
            continue

        # Dead cell decimation
        if type_id == 7 and DEAD_DECIMATE != 1:
            if DEAD_DECIMATE == 0:
                print(f"  {spec['name']}: {count:,} particles -- SKIPPED (DEAD_DECIMATE=0)")
                continue

            indices = np.where(mask)[0]
            # Keep every Nth dead cell
            indices = indices[::DEAD_DECIMATE]
            mask = np.zeros(len(types), dtype=bool)
            mask[indices] = True
            kept = mask.sum()
            print(f"  {spec['name']}: {count:,} -> {kept:,} particles "
                  f"(decimated 1/{DEAD_DECIMATE})")
            count = kept
        else:
            print(f"  {spec['name']}: {count:,} particles")

        # Extract this species' data
        pos = positions_scaled[mask]
        rad = radii_scaled[mask]

        # Create point cloud object
        obj_name = f"Biofilm_{spec['name']}"
        obj = create_point_cloud(obj_name, pos, rad)
        total_imported += count

        # Create material
        material = create_material(type_id)

        # Set up Geometry Nodes
        setup_geometry_nodes(obj, sphere_obj, material)

        created_objects.append(obj)

    print(f"\nTotal imported: {total_imported:,} particles")
    print(f"Objects created: {len(created_objects)}")

    # --- Scene setup ---
    if SETUP_SCENE:
        print("\nSetting up scene...")
        setup_scene(domain_center, domain_size)

    # --- Organize into collection ---
    collection = bpy.data.collections.new("Mercury_Reef")
    bpy.context.scene.collection.children.link(collection)

    for obj in created_objects:
        # Move from default collection to Mercury_Reef
        for coll in obj.users_collection:
            coll.objects.unlink(obj)
        collection.objects.link(obj)

    # Move sphere template too
    for coll in sphere_obj.users_collection:
        coll.objects.unlink(sphere_obj)
    collection.objects.link(sphere_obj)

    print("\n" + "=" * 60)
    print("  IMPORT COMPLETE")
    print("=" * 60)
    print(f"""
Next steps:
  1. In the viewport, press Z -> Material Preview to see colors
  2. Each species is a separate object -- toggle visibility in the Outliner
  3. DEAD cells are the heaviest layer -- hide them first for fast preview
  4. To render: F12 (or Render -> Render Image)
  5. For production quality:
     - Increase Cycles samples to 1024-2048
     - Add an HDRI to World shader (download from polyhaven.com)
     - Adjust camera position and depth of field to taste
     - Increase SPHERE_SUBDIVISIONS to 3

Estimated render time:
  - Preview (256 samples, 1080p): ~2-5 min (GPU) / ~10-30 min (CPU)
  - Production (1024 samples, 4K): ~15-60 min (GPU) / ~2-8 hrs (CPU)
""")


# ╔══════════════════════════════════════════════════════════════╗
# ║                   ANIMATION LOADER                          ║
# ║  (Optional — uncomment and run separately for time series)  ║
# ╚══════════════════════════════════════════════════════════════╝

def load_animation(vtk_dir, step_interval=80, frame_rate=15):
    """
    Load a time series of VTU files as a Blender animation.

    This replaces the mesh data on each frame using a frame_change handler.
    Call this AFTER running main() to set up materials and Geometry Nodes.

    Usage:
        load_animation(r"C:\\Users\\ltjjp\\OneDrive\\Desktop\\mercury_reef_vtk_combined")

    Args:
        vtk_dir: Directory containing dump*.vtu files
        step_interval: NUFEB steps between VTU dumps (default 80)
        frame_rate: Blender playback frame rate
    """
    import glob
    import re

    # Find all VTU files and sort numerically
    vtu_files = glob.glob(os.path.join(vtk_dir, "dump*.vtu"))
    vtu_files.sort(key=lambda f: int(re.search(r'dump(\d+)', f).group(1)))

    print(f"Found {len(vtu_files)} VTU files for animation")

    # Store file list as a scene property
    bpy.context.scene["vtu_files"] = [str(f) for f in vtu_files]
    bpy.context.scene["vtu_scale"] = SCALE
    bpy.context.scene["vtu_dead_decimate"] = DEAD_DECIMATE

    # Set timeline
    bpy.context.scene.frame_start = 0
    bpy.context.scene.frame_end = len(vtu_files) - 1
    bpy.context.scene.render.fps = frame_rate

    def update_frame(scene):
        """Frame change handler: reload VTU data for current frame."""
        frame = scene.frame_current
        files = scene.get("vtu_files", [])

        if frame < 0 or frame >= len(files):
            return

        filepath = files[frame]
        scale = scene.get("vtu_scale", 1e6)
        dead_dec = scene.get("vtu_dead_decimate", 5)

        try:
            positions, types_arr, diameters = parse_vtu(filepath)
        except Exception as e:
            print(f"Error loading frame {frame}: {e}")
            return

        positions_scaled = positions * scale
        radii_scaled = (diameters * scale) / 2.0

        # Update each species object
        for type_id, spec in SPECIES.items():
            obj_name = f"Biofilm_{spec['name']}"
            obj = bpy.data.objects.get(obj_name)
            if obj is None:
                continue

            mask = types_arr == type_id
            if type_id == 7 and dead_dec not in (1, None):
                if dead_dec == 0:
                    mask[:] = False
                else:
                    indices = np.where(mask)[0][::dead_dec]
                    mask = np.zeros(len(types_arr), dtype=bool)
                    mask[indices] = True

            pos = positions_scaled[mask]
            rad = radii_scaled[mask]

            # Rebuild mesh
            mesh = obj.data
            mesh.clear_geometry()

            if len(pos) > 0:
                mesh.vertices.add(len(pos))
                mesh.vertices.foreach_set("co", pos.flatten().tolist())

                # Update radius attribute
                if "radius" in mesh.attributes:
                    mesh.attributes.remove(mesh.attributes["radius"])
                attr = mesh.attributes.new(
                    name="radius", type='FLOAT', domain='POINT'
                )
                attr.data.foreach_set("value", rad.tolist())

            mesh.update()

    # Register the handler
    # Remove any existing handler first
    for handler in bpy.app.handlers.frame_change_post:
        if handler.__name__ == "update_frame":
            bpy.app.handlers.frame_change_post.remove(handler)

    bpy.app.handlers.frame_change_post.append(update_frame)
    print(f"Animation handler registered: {len(vtu_files)} frames")
    print("Use the timeline to scrub through the simulation.")
    print("NOTE: Each frame reloads a VTU file -- scrubbing will have a delay.")


# --- Run ---
if __name__ == "__main__":
    main()

    # Uncomment the line below to enable animation (time series):
    load_animation(r"C:\Users\ltjjp\OneDrive\Desktop\mercury_reef_vtk_combined")
