"""
Script: Visualize simplified PRg structure using PyVista
--------------------------------------------------------
This script renders a simplified 3D model of a polyrotaxane glass (PRg)
using PyVista. It displays PEG chains (as tubes) and CD molecules (as hollow
frustums) oriented along extracted normal vectors.

Input files:
    1. Backbone atoms file          (e.g., 2_backbone_xxx.data)
    2. CD centers coordinate file   (e.g., 3_cd_centers_xxx.data)
    3. CD normal vectors file       (e.g., 4_cd_normals_xxx.data)

Output:
    - PNG image (for paper figures or animation)
    - (Optional) Interactive 3D view if enabled

Usage:
    python3 5_plot.py backbone.data cd_centers.data cd_normals.data output.png

Notes:
    - To enable interactive viewing:
        → Set USE_INTERACTIVE_VIEW = True
        → Change: off_screen=False
        → Uncomment: plotter.show(auto_close=False)
    - GPU acceleration (optional):
        export __NV_PRIME_RENDER_OFFLOAD=1
        export __GLX_VENDOR_LIBRARY_NAME=nvidia

"""

import os
import sys
import numpy as np
import pyvista as pv


# ====================== User Settings ======================
USE_INTERACTIVE_VIEW = False  # True → open interactive window

config = {
    # --- Cyclodextrin (CD) settings ---
    'r1_outer': 7.75,         # outer radius (bottom)
    'r2_outer': 6.0,          # outer radius (top)
    'r_inner': 4.0,           # inner radius
    'h': 13.0,                # height of truncated cone
    'num_sides': 90,          # circular subdivision
    'transparency': 0.5,
    'CD_color': '#f2d236',

    # --- Chain and end-cap ---
    'tube_radius': 1.5,
    'chain_color': '#5b63f3',
    'sphere_radius': 4.0,
    'sphere_color': '#bb2638',
    'sphere_resolution': 50,
    'num_sides_tube': 40,
    'distance_threshold': 50,   # skip broken links
    'show_individual_ends': True,
    'use_individual_colors': False,

    # --- Lighting and material ---
    'ambient': 0.2,
    'diffuse': 0.8,
    'specular': 0.1,
    'roughness': 1.0,
    'ambient_CD': 0.4,
    'diffuse_CD': 0.8,
    'specular_CD': 0.5,
    'roughness_CD': 1.0,

    # --- Rendering window ---
    'width': 3840,
    'height': 2160,
    'x_span': 0,
    'y_span': 0,
    'z_span': 1,
}
# =============================================================


# ====================== Data Import Functions ======================
def import_backbone_data(file_path):
    """
    Read backbone chain coordinates from 'molecule_x' formatted file.
    Returns:
        chains: list of chains (list of coordinates)
        chain_ends: list of start & end points
    """
    chains, chain_ends, all_chain_ends = [], [], []
    current_chain = []

    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("molecule_"):
                if current_chain:
                    chains.append(current_chain)
                    chain_ends.append({
                        'atoms': [current_chain[0], current_chain[-1]]
                    })
                    all_chain_ends.extend([current_chain[0], current_chain[-1]])
                    current_chain = []
            else:
                parts = line.split()
                if len(parts) >= 6:
                    x, y, z = map(float, parts[3:6])
                    current_chain.append((x, y, z))

    if current_chain:
        chains.append(current_chain)
        chain_ends.append({'atoms': [current_chain[0], current_chain[-1]]})
        all_chain_ends.extend([current_chain[0], current_chain[-1]])

    return chains, chain_ends, all_chain_ends


def load_cd_centers(file_path):
    """Load CD molecule centers as {mol_id: (x, y, z)}."""
    centers = {}
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) == 4:
                mol_id, x, y, z = int(parts[0]), *map(float, parts[1:])
                centers[mol_id] = (x, y, z)
    return centers


def load_cd_normals(file_path):
    """Load CD normal vectors as {mol_id: (vx, vy, vz)}."""
    normals = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.split()
            if len(parts) == 4:
                mol_id, vx, vy, vz = int(parts[0]), *map(float, parts[1:])
                normals[mol_id] = (vx, vy, vz)
    return normals


# ====================== Math Utilities ======================
def rotation_matrix_from_vectors(vec1, vec2):
    """Return rotation matrix that rotates vec1 to vec2."""
    a = vec1 / np.linalg.norm(vec1)
    b = vec2 / np.linalg.norm(vec2)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)

    if np.isclose(s, 0):
        if c > 0:  # parallel
            return np.eye(3)
        # anti-parallel
        axis = np.array([1, 0, 0]) if not np.allclose(a[:2], [0, 0]) else np.array([0, 1, 0])
        v = np.cross(a, axis) / np.linalg.norm(np.cross(a, axis))
        kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        return np.eye(3) + 2 * kmat.dot(kmat)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))


# ====================== Geometry Construction ======================
def create_hollow_frustum(r1_outer, r2_outer, r_inner, h, num_sides=90):
    """Generate a hollow truncated cone mesh (CD ring)."""
    theta = np.linspace(0, 2 * np.pi, num_sides, endpoint=False)
    ring0 = np.column_stack([r1_outer * np.cos(theta), r1_outer * np.sin(theta), np.zeros(num_sides)])
    ring1 = np.column_stack([r2_outer * np.cos(theta), r2_outer * np.sin(theta), np.full(num_sides, h)])
    ring2 = np.column_stack([r_inner * np.cos(theta), r_inner * np.sin(theta), np.full(num_sides, h)])
    ring3 = np.column_stack([r_inner * np.cos(theta), r_inner * np.sin(theta), np.zeros(num_sides)])

    vertices = np.vstack([ring0, ring1, ring2, ring3])
    n = num_sides
    faces = []
    for i in range(n):
        ni = (i + 1) % n
        faces.extend([
            [3, i, ni, n + i],
            [3, ni, n + ni, n + i],
            [3, n + i, n + ni, 2 * n + i],
            [3, n + ni, 2 * n + ni, 2 * n + i],
            [3, 3 * n + i, 3 * n + ni, 2 * n + i],
            [3, 3 * n + ni, 2 * n + ni, 2 * n + i],
            [3, i, ni, 3 * n + i],
            [3, ni, 3 * n + ni, 3 * n + i],
        ])
    return vertices, np.array(faces).flatten().astype(np.int32)


def create_tube_from_chain(chain_coords, radius=1.0, num_sides=20, distance_threshold=50):
    """Generate a tube mesh along a polymer chain."""
    chain = np.array(chain_coords)
    theta = np.linspace(0, 2 * np.pi, num_sides, endpoint=False)
    circle_x, circle_y = np.cos(theta), np.sin(theta)
    verts = []

    for i, atom in enumerate(chain):
        if i == 0:
            direction = chain[i + 1] - atom
        elif i == len(chain) - 1:
            direction = atom - chain[i - 1]
        else:
            direction = (chain[i + 1] - atom) + (atom - chain[i - 1])
        direction /= np.linalg.norm(direction)

        orth1 = np.cross(direction, [0, 0, 1]) if not np.allclose(direction, [0, 0, 1]) else np.array([1, 0, 0])
        orth1 /= np.linalg.norm(orth1)
        orth2 = np.cross(direction, orth1)
        for t in range(num_sides):
            verts.append(atom + radius * (circle_x[t] * orth1 + circle_y[t] * orth2))

    verts = np.array(verts)
    faces = []
    for seg in range(len(chain) - 1):
        if np.linalg.norm(chain[seg + 1] - chain[seg]) > distance_threshold:
            continue
        b, n = seg * num_sides, (seg + 1) * num_sides
        for t in range(num_sides):
            tn = (t + 1) % num_sides
            faces.extend([[3, b + t, b + tn, n + t], [3, b + tn, n + tn, n + t]])
    return verts, np.array(faces).flatten().astype(np.int32)


# ====================== Main Visualization ======================
def main():
    input_chain, input_center, input_normal, output_png = sys.argv[1:5]
    current_dir = os.path.dirname(os.path.abspath(__file__))
    chain_file = os.path.join(current_dir, input_chain)
    center_file = os.path.join(current_dir, input_center)
    normal_file = os.path.join(current_dir, input_normal)

    chains, chain_ends, _ = import_backbone_data(chain_file)
    centers = load_cd_centers(center_file)
    normals = load_cd_normals(normal_file)

    plotter = pv.Plotter(
        window_size=(config['width'], config['height']),
        off_screen=not USE_INTERACTIVE_VIEW
    )
    plotter.background_color = "white"

    # --- (A) Draw PEG chains ---
    for i, chain in enumerate(chains):
        verts, faces = create_tube_from_chain(
            chain, radius=config['tube_radius'],
            num_sides=config['num_sides_tube'],
            distance_threshold=config['distance_threshold']
        )
        tube = pv.PolyData(verts, faces)
        color = np.random.rand(3) if config['use_individual_colors'] else config['chain_color']
        plotter.add_mesh(tube, color=color, opacity=1.0, smooth_shading=True,
                         specular=config['specular'], diffuse=config['diffuse'],
                         ambient=config['ambient'], roughness=config['roughness'])

        if config['show_individual_ends']:
            for end in chain_ends[i]['atoms']:
                sphere = pv.Sphere(radius=config['sphere_radius'],
                                   center=end,
                                   theta_resolution=config['sphere_resolution'],
                                   phi_resolution=config['sphere_resolution'])
                plotter.add_mesh(sphere, color=config['sphere_color'],
                                 specular=config['specular'], diffuse=config['diffuse'],
                                 ambient=config['ambient'], roughness=config['roughness'])

    # --- (B) Draw CDs ---
    for mol_id, center in centers.items():
        vx, vy, vz = normals[mol_id]
        verts, faces = create_hollow_frustum(
            config['r1_outer'], config['r2_outer'], config['r_inner'],
            config['h'], config['num_sides']
        )
        cd_mesh = pv.PolyData(verts, faces)
        R = rotation_matrix_from_vectors(np.array([0, 0, 1]), np.array([vx, vy, vz]))
        cd_mesh.points = np.dot(cd_mesh.points, R.T)
        shift = np.array(center) - cd_mesh.points.mean(axis=0)
        cd_mesh.points += shift

        plotter.add_mesh(cd_mesh, color=config['CD_color'],
                         opacity=config['transparency'], smooth_shading=True,
                         specular=config['specular_CD'], diffuse=config['diffuse_CD'],
                         ambient=config['ambient_CD'], roughness=config['roughness_CD'])

    # --- (C) Camera and rendering ---
    try:
        plotter.enable_depth_peeling()
    except:
        pass
    plotter.enable_parallel_projection()
    plotter.camera_position = [
        (config['x_span'], config['y_span'], config['z_span']),
        (0, 0, 0),
        (0, 1, 0)
    ]
    plotter.reset_camera()
    plotter.camera.parallel_scale = 200

    # --- (D) Output ---
    if USE_INTERACTIVE_VIEW:
        plotter.show(auto_close=False)
    else:
        output_path = os.path.join(current_dir, output_png)
        plotter.screenshot(output_path)
        print(f"✅ PNG image saved to: {output_path}")
        plotter.close()


if __name__ == "__main__":
    main()
