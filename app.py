# app.py
import os
import sys
import glob
import numpy as np
import trimesh
import plotly.graph_objs as go
import gradio as gr
from plyfile import PlyData

# ------------------ Load pyregister module ------------------
so_candidates = glob.glob("build/python/pyregister*.so")
if not so_candidates:
    os.system("bash build.sh")
    so_candidates = glob.glob("build/python/pyregister*.so")
if not so_candidates:
    raise FileNotFoundError("pyregister build failed, no .so found in build/python")
sys.path.append("build/python")

import pyregister  

# ------------------ Example data ------------------
EXAMPLES = {
    "Rigid":    {"target": "./examples/data/fricp/target.ply", "source": "./examples/data/fricp/source.ply"},
    "NonRigid": {"target": "./examples/data/spare/target.obj", "source": "./examples/data/spare/source.obj"},
}

FRICP_PARAS = {
    "useinit": False,
    "fileinit": "",
    "maxiter": 100,
    "stop": 1e-5
}

SPARE_PARAS = {
    "iters": 30,
    "stopcoarse": 1e-3,
    "stopfine": 1e-4,
    "use_landmark": False,
    "src": [],
    "tar": []
}

def load_example(example_type):
    target_path = EXAMPLES[example_type]["target"]
    source_path = EXAMPLES[example_type]["source"]
    method = example_type
    return source_path, target_path, method

def read_mesh(file_path):
    ext = os.path.splitext(file_path)[1].lower()
    if ext == ".ply":
        plydata = PlyData.read(file_path)
        vertex = plydata["vertex"]
        vertices = np.stack([vertex["x"], vertex["y"], vertex["z"]], axis=-1)
        faces = None
        if "face" in plydata and len(plydata["face"].data) > 0:
            face_list = [f[0] for f in plydata["face"].data]
            if len(face_list) > 0:
                faces = np.array(face_list)
                if faces.ndim == 1:
                    faces = faces.reshape(-1, 3)
        return vertices, faces
    elif ext == ".obj":
        mesh = trimesh.load(file_path, process=False)
        return np.asarray(mesh.vertices), (np.asarray(mesh.faces) if mesh.faces.size > 0 else None)
    else:
        raise ValueError(f"Unsupported file type: {ext}")

def plot_source_target_mesh(target_file, source_file, alpha_target=1.0, alpha_source=1.0, scatter_mode=False):
    target_v, target_f = read_mesh(target_file)
    source_v, source_f = read_mesh(source_file)

    fig = go.Figure()
    lighting_opts = dict(ambient=0.5, diffuse=0.8, specular=0.6, roughness=0.25)
    light_pos = dict(x=100, y=200, z=50)

    target_color = "crimson"

    if scatter_mode:

        fig.add_trace(go.Scatter3d(
            x=target_v[:, 0], y=target_v[:, 1], z=target_v[:, 2],
            mode="markers",
  
            marker=dict(size=0.6, color=target_color, opacity=1.0),
            name="Target Points",
            hoverinfo="text",
            text=[
                f"<b><span style='color:black;'>Target ID:</span></b> {i}"
                f"<br><span style='font-size:12px;'>x={x:.3f}, y={y:.3f}, z={z:.3f}</span>"
                for i, (x, y, z) in enumerate(target_v)
            ],
            hoverlabel=dict(
                bgcolor="white",
                bordercolor="black",
                font=dict(color="black", size=13, family="Arial", weight="bold")
            )
        ))

        fig.add_trace(go.Scatter3d(
            x=source_v[:, 0], y=source_v[:, 1], z=source_v[:, 2],
            mode="markers",

            marker=dict(size=0.6, color="limegreen", opacity=1.0),
            name="Source Points",
            hoverinfo="text",
            text=[
                f"<b><span style='color:red;'>Source ID:</span></b> {i}"
                f"<br><span style='font-size:12px;'>x={x:.3f}, y={y:.3f}, z={z:.3f}</span>"
                for i, (x, y, z) in enumerate(source_v)
            ],
            hoverlabel=dict(
                bgcolor="white",
                bordercolor="limegreen",
                font=dict(color="red", size=13, family="Arial", weight="bold")
            )
        ))

    else:

        if target_f is not None:
            fig.add_trace(go.Mesh3d(
                x=target_v[:, 0], y=target_v[:, 1], z=target_v[:, 2],
                i=target_f[:, 0], j=target_f[:, 1], k=target_f[:, 2],
                color="khaki", opacity=alpha_target, name="Target Mesh",
                lighting=lighting_opts, lightposition=light_pos,
                hovertemplate="<b>Target Surface</b><br>x=%{x:.3f}<br>y=%{y:.3f}<br>z=%{z:.3f}<extra></extra>",
            ))
        else:
            fig.add_trace(go.Scatter3d(
                x=target_v[:, 0], y=target_v[:, 1], z=target_v[:, 2],
                mode="markers",
                marker=dict(size=0.6, color=target_color, opacity=alpha_target),  
                name="Target Points"
            ))

        if source_f is not None:
            fig.add_trace(go.Mesh3d(
                x=source_v[:, 0], y=source_v[:, 1], z=source_v[:, 2],
                i=source_f[:, 0], j=source_f[:, 1], k=source_f[:, 2],
                color="darkseagreen", opacity=alpha_source, name="Source Mesh",
                lighting=lighting_opts, lightposition=light_pos,
                hovertemplate="<b>Source Surface</b><br>x=%{x:.3f}<br>y=%{y:.3f}<br>z=%{z:.3f}<extra></extra>",
            ))
        else:
            fig.add_trace(go.Scatter3d(
                x=source_v[:, 0], y=source_v[:, 1], z=source_v[:, 2],
                mode="markers",
                marker=dict(size=0.6, color="limegreen", opacity=alpha_source),  
                name="Source Points"
            ))

    fig.update_layout(
        height=600, width=600, margin=dict(l=0, r=0, t=40, b=0),
        scene=dict(aspectmode="data"), showlegend=True
    )
    return fig



def plot_result_target_mesh(target_file, result_file):
    target_v, target_f = read_mesh(target_file)
    result_v, result_f = read_mesh(result_file)

    fig = go.Figure()
    lighting_opts = dict(ambient=0.5, diffuse=0.8, specular=0.6, roughness=0.25)
    light_pos = dict(x=100, y=200, z=50)

    target_mesh_color = "khaki"
    result_mesh_color = "lightblue"

    if target_f is not None:
        fig.add_trace(go.Mesh3d(
            x=target_v[:, 0], y=target_v[:, 1], z=target_v[:, 2],
            i=target_f[:, 0], j=target_f[:, 1], k=target_f[:, 2],
            color=target_mesh_color, opacity=1.0, name="Target Mesh",
            lighting=lighting_opts, lightposition=light_pos
        ))
    else:
        fig.add_trace(go.Scatter3d(
            x=target_v[:, 0], y=target_v[:, 1], z=target_v[:, 2],
            mode="markers",
            marker=dict(size=0.6, color="crimson", opacity=1.0),
            name="Target Points (Red)"
        ))

    if result_f is not None:
        fig.add_trace(go.Mesh3d(
            x=result_v[:, 0], y=result_v[:, 1], z=result_v[:, 2],
            i=result_f[:, 0], j=result_f[:, 1], k=result_f[:, 2],
            color=result_mesh_color, opacity=1.0, name="Result Mesh",
            lighting=lighting_opts, lightposition=light_pos
        ))
    else:
        fig.add_trace(go.Scatter3d(
            x=result_v[:, 0], y=result_v[:, 1], z=result_v[:, 2],
            mode="markers",
            marker=dict(size=0.6, color="royalblue", opacity=1.0),
            name="Result Points (Blue)"
        ))

    fig.update_layout(height=600, width=600, margin=dict(l=0, r=0, t=40, b=0),
                      scene=dict(aspectmode="data"), showlegend=True)
    return fig

# ------------------ Registration ------------------
def register_and_visualize_with_zip(target_file, source_file, output_dir, method,
                                    useinit=False, fileinit=None, maxiter=100, stop=1e-5,
                                    iters=30, stopcoarse=1e-3, stopfine=1e-4):
    if target_file is None or source_file is None:
        raise gr.Error("Please upload both target and source point cloud files first!")

    target_path = target_file.name
    source_path = source_file.name
    os.makedirs(output_dir, exist_ok=True)

    # === Rigid ===
    if method == "Rigid":
        reg = pyregister.RigidFricpRegistration()
        FRICP_PARAS["useinit"] = bool(useinit)
        FRICP_PARAS["fileinit"] = fileinit.name if useinit and fileinit else ""
        FRICP_PARAS["maxiter"] = int(maxiter)
        FRICP_PARAS["stop"] = float(stop)

        reg.Paras_init(FRICP_PARAS["useinit"],
                       FRICP_PARAS["fileinit"],
                       FRICP_PARAS["maxiter"],
                       FRICP_PARAS["stop"])
        output_file = "FRICP_res.ply"

    # === NonRigid ===
    else:
        reg = pyregister.NonrigidSpareRegistration()
        SPARE_PARAS["iters"] = int(iters)
        SPARE_PARAS["stopcoarse"] = float(stopcoarse)
        SPARE_PARAS["stopfine"] = float(stopfine)

        reg.Paras_init(
            SPARE_PARAS["iters"],
            SPARE_PARAS["stopcoarse"],
            SPARE_PARAS["stopfine"],
            SPARE_PARAS["use_landmark"],
            SPARE_PARAS["src"],
            SPARE_PARAS["tar"]
        )
        output_file = "spare_res.ply"


    reg.Reg(target_path, source_path, output_dir)

    result_path = os.path.join(output_dir, output_file)
    fig_result = plot_result_target_mesh(target_path, result_path)

    import zipfile
    zip_path = os.path.join(output_dir, "results.zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
        for f in os.listdir(output_dir):
            if f.endswith(".ply"):
                zipf.write(os.path.join(output_dir, f), arcname=f)

    return fig_result, zip_path


# ------------------ Reregister ------------------
def reuse_last_result_as_source(target_file, output_dir):
    if target_file is None:
        raise gr.Error("Please upload the target point cloud file first!")
    if not os.path.exists(output_dir):
        raise gr.Error("Output directory does not exist!")

    result_candidates = []
    for fname in ["FRICP_res.ply", "spare_res.ply"]:
        fpath = os.path.join(output_dir, fname)
        if os.path.exists(fpath):
            result_candidates.append(fpath)
    if not result_candidates:
        raise gr.Error("No previous registration result found. Please run registration first!")

    result_candidates.sort(key=os.path.getmtime, reverse=True)
    latest_result = result_candidates[0]

    fig = plot_source_target_mesh(target_file.name, latest_result)
    return latest_result, fig

# ------------------ Utility ------------------
def reset_parameters():
    """
    Reset all registration parameters (Rigid + NonRigid + Landmarks)
    and clear any manual inputs in the UI.
    This also resets the 'Use Landmarks' checkbox and hides the landmark UI group.
    """

    # --- Reset Rigid (FRICP) parameters ---
    FRICP_PARAS.update({
        "useinit": False,
        "fileinit": "",
        "maxiter": 100,
        "stop": 1e-5
    })

    # --- Reset NonRigid (Spare) parameters ---
    SPARE_PARAS.update({
        "iters": 30,
        "stopcoarse": 1e-3,
        "stopfine": 1e-4,
        "use_landmark": False,
        "src": [],
        "tar": []
    })

    print("[reset] All parameters reset.")
    print("[reset] FRICP_PARAS:", FRICP_PARAS)
    print("[reset] SPARE_PARAS:", SPARE_PARAS)


    return (
        FRICP_PARAS["useinit"],     # Checkbox
        None,                       # File input reset
        FRICP_PARAS["maxiter"],
        FRICP_PARAS["stop"],
        SPARE_PARAS["iters"],
        SPARE_PARAS["stopcoarse"],
        SPARE_PARAS["stopfine"],
        gr.update(value=""),        
        gr.update(value=""),      
        gr.update(value=False),    
        gr.update(visible=False)   
    )


def clear_all():
    """
    Full reset for all UI elements:
    - Clears uploaded files
    - Resets plots
    - Resets all parameters (Rigid + NonRigid + Landmark)
    - Resets landmark checkbox & hides group
    """

    FRICP_PARAS.update({
        "useinit": False,
        "fileinit": "",
        "maxiter": 100,
        "stop": 1e-5
    })
    SPARE_PARAS.update({
        "iters": 30,
        "stopcoarse": 1e-3,
        "stopfine": 1e-4,
        "use_landmark": False,
        "src": [],
        "tar": []
    })

    print("[clear] All files, parameters, and landmarks cleared.")
    return (
        None,  # target_input
        None,  # source_input
        None,  # upload_plot
        None,  # result_plot
        FRICP_PARAS["useinit"],  # fricp_useinit
        None,  # fricp_fileinit
        FRICP_PARAS["maxiter"],
        FRICP_PARAS["stop"],
        SPARE_PARAS["iters"],
        SPARE_PARAS["stopcoarse"],
        SPARE_PARAS["stopfine"],
        gr.update(value=""),    
        gr.update(value=""),      
        gr.update(value=False),   
        gr.update(visible=False)  
    )


def visualize_and_store(target_file, source_file, scatter_mode=False):
    if target_file is None or source_file is None:
        return None, None, None
    target_v, _ = read_mesh(target_file.name)
    source_v, _ = read_mesh(source_file.name)
    fig = plot_source_target_mesh(target_file.name, source_file.name, scatter_mode=scatter_mode)

    return fig, source_v, target_v

def clear_landmarks():
    return [], [], "", ""
def highlight_landmarks_on_mesh(target_file, source_file, src_text, tar_text):
    src_v, _ = read_mesh(source_file.name)
    tar_v, _ = read_mesh(target_file.name)
    src_ids = [int(i) for i in src_text.split(",") if i.strip().isdigit()]
    tar_ids = [int(i) for i in tar_text.split(",") if i.strip().isdigit()]
   

    fig = plot_source_target_mesh(target_file.name, source_file.name)

    # === Source landmarks ===
    if len(src_ids) > 0:
        pts = src_v[src_ids]
        fig.add_trace(go.Scatter3d(
            x=pts[:, 0], y=pts[:, 1], z=pts[:, 2],
            mode="markers+text",
            text=[str(i) for i in src_ids],
            textposition="top center",
            textfont=dict(color="black", size=10, family="Arial"),
            marker=dict(
                size=6,
                color="limegreen",
                opacity=0.9,
                line=dict(width=2, color="white"),
                symbol="circle"
            ),
            name="Source Landmarks",
            hoverinfo="text",
            hovertext=[
                f"<b>Source ID:</b> {i}<br>x={x:.3f}, y={y:.3f}, z={z:.3f}"
                for i, (x, y, z) in zip(src_ids, pts)
            ]
        ))

    # === Target landmarks ===
    if len(tar_ids) > 0:
        pts = tar_v[tar_ids]
        fig.add_trace(go.Scatter3d(
            x=pts[:, 0], y=pts[:, 1], z=pts[:, 2],
            mode="markers+text",
            text=[str(i) for i in tar_ids],
            textposition="top center",
            textfont=dict(color="black", size=10, family="Arial"),
            marker=dict(
                size=6,
                color="crimson",
                opacity=0.9,
                line=dict(width=2, color="white"),
                symbol="circle"
            ),
            name="Target Landmarks",
            hoverinfo="text",
            hovertext=[
                f"<b>Target ID:</b> {i}<br>x={x:.3f}, y={y:.3f}, z={z:.3f}"
                for i, (x, y, z) in zip(tar_ids, pts)
            ]
        ))


    fig.update_layout(
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False)
        ),
        legend=dict(bgcolor="rgba(255,255,255,0.6)"),
        margin=dict(l=0, r=0, t=40, b=0),
        title="Landmark Highlight View"
    )

    return fig

def start_landmark_selection(target_file, source_file):
    if target_file is None or source_file is None:
        raise gr.Error("Please upload both Source and Target first!")

    fig = plot_source_target_mesh(target_file.name, source_file.name, scatter_mode=True)
    return fig
def update_landmark_ids(src_text, tar_text):
    src_ids = [int(i) for i in src_text.split(",") if i.strip().isdigit()]
    tar_ids = [int(i) for i in tar_text.split(",") if i.strip().isdigit()]
    return src_ids, tar_ids, src_text, tar_text


def set_landmarks_to_registration(src_ids, tar_ids):
    if not src_ids or not tar_ids:
        raise gr.Error("Please select both Source and Target points first!")
    if len(src_ids) != len(tar_ids):
        raise gr.Error(f"Landmark count mismatch: Source={len(src_ids)} vs Target={len(tar_ids)}")

    SPARE_PARAS["use_landmark"] = True
    SPARE_PARAS["src"] = [int(i) for i in src_ids]
    SPARE_PARAS["tar"] = [int(i) for i in tar_ids]

    return f"✅ Landmarks staged: {len(src_ids)} Source ↔ {len(tar_ids)} Target. Will be used in NonRigid registration."

# ------------------ UI ------------------
CSS_FIX = """
html, body, .gradio-container {
    min-height: 100vh;
    overflow-y: auto;
}
"""



with gr.Blocks(css=CSS_FIX) as demo:
    gr.Markdown("## Point Cloud / Mesh Registration Visualization Tool")

    src_ids_state = gr.State([])
    tar_ids_state = gr.State([])
    src_vertices_state = gr.State(None)
    tar_vertices_state = gr.State(None)

    with gr.Row():
        with gr.Column(scale=1):

            source_input = gr.File(label="Source File", file_types=[".ply", ".obj"])
            target_input = gr.File(label="Target File", file_types=[".ply", ".obj"])
            method_dropdown = gr.Dropdown(label="Registration Method",
                                          choices=["Rigid", "NonRigid"], value="Rigid")

            with gr.Accordion("Rigid Parameters", open=False):
                fricp_useinit = gr.Checkbox(label="Use initial transform?", value=False)
                fricp_fileinit = gr.File(label="Initial transform file", file_types=[".txt"])
                fricp_maxiter = gr.Number(label="Max iterations", value=100)
                fricp_stop = gr.Number(label="Stop threshold", value=1e-5)


            with gr.Accordion("NonRigid Parameters", open=False):
                spare_iters = gr.Number(label="Iterations", value=30)
                spare_stopcoarse = gr.Number(label="Stop Coarse", value=1e-3)
                spare_stopfine = gr.Number(label="Stop Fine", value=1e-4)


                use_landmark_cb = gr.Checkbox(label="Use Landmarks", value=False)


                with gr.Group(visible=False) as landmark_group:
                    gr.Markdown("### Landmark Point Selection")
                    selected_src = gr.Textbox(
                        label="Selected Source IDs (comma-separated)",
                        placeholder="e.g. 12,45,89", interactive=True)
                    selected_tar = gr.Textbox(
                        label="Selected Target IDs (comma-separated)",
                        placeholder="e.g. 7,42,105", interactive=True)

                    with gr.Row():
                        start_landmark_btn = gr.Button("Start Landmark Selection", elem_classes="primary")
                        highlight_btn = gr.Button("Highlight Landmarks")
                        clear_landmark_btn = gr.Button("🧹 Clear Selections")

                    select_button = gr.Button("Confirm Landmarks to Registration", elem_classes="primary")


            output_dir = gr.Textbox(label="Output Directory", value="./output/", placeholder="./output/")
            with gr.Row():
                run_button = gr.Button("Run Registration", elem_classes="primary")
                rerun_button = gr.Button("Reregister (Use Last Result as Source)", elem_classes="primary")
                clear_button = gr.Button("Clear", elem_classes="secondary")
                reset_button = gr.Button("Reset Parameters", elem_classes="primary")

            
            with gr.Row():
                example_dropdown = gr.Dropdown(label="Load Example Data",
                                               choices=["Rigid", "NonRigid"], value="Rigid")
                example_button = gr.Button("Load Example")

            example_button.click(fn=load_example,
                                 inputs=[example_dropdown],
                                 outputs=[source_input, target_input, method_dropdown])

        
        with gr.Column(scale=2):
            upload_plot = gr.Plot(label="Source vs Target Mesh")
            result_plot = gr.Plot(label="Result vs Target Mesh")
            download_button = gr.File(label="Download Result Directory",
                                      file_types=[".zip"], interactive=False)


    reset_button.click(
        fn=reset_parameters,
        inputs=None,
        outputs=[
            fricp_useinit, fricp_fileinit, fricp_maxiter, fricp_stop,
            spare_iters, spare_stopcoarse, spare_stopfine,
            selected_src, selected_tar,
            use_landmark_cb, landmark_group
        ]
    )

    def toggle_landmark_visibility(use_landmark):
        SPARE_PARAS["use_landmark"] = bool(use_landmark)
        print(f"[UI] use_landmark set to {SPARE_PARAS['use_landmark']}")
        return gr.update(visible=use_landmark)

    use_landmark_cb.change(fn=toggle_landmark_visibility,
                           inputs=[use_landmark_cb],
                           outputs=[landmark_group])

    clear_button.click(
    fn=clear_all,
    inputs=None,
    outputs=[
        target_input, source_input, upload_plot, result_plot,
        fricp_useinit, fricp_fileinit, fricp_maxiter, fricp_stop,
        spare_iters, spare_stopcoarse, spare_stopfine,
        selected_src, selected_tar,
        use_landmark_cb, landmark_group
    ]
)


    target_input.change(fn=visualize_and_store,
                        inputs=[target_input, source_input],
                        outputs=[upload_plot, src_vertices_state, tar_vertices_state])
    source_input.change(fn=visualize_and_store,
                        inputs=[target_input, source_input],
                        outputs=[upload_plot, src_vertices_state, tar_vertices_state])

    start_landmark_btn.click(fn=start_landmark_selection,
                             inputs=[target_input, source_input],
                             outputs=[upload_plot])

    highlight_btn.click(fn=highlight_landmarks_on_mesh,
                        inputs=[target_input, source_input, selected_src, selected_tar],
                        outputs=[upload_plot])

    clear_landmark_btn.click(fn=clear_landmarks,
                             inputs=None,
                             outputs=[src_ids_state, tar_ids_state, selected_src, selected_tar])

    select_button.click(fn=update_landmark_ids,
                        inputs=[selected_src, selected_tar],
                        outputs=[src_ids_state, tar_ids_state, selected_src, selected_tar]
                        ).then(fn=set_landmarks_to_registration,
                               inputs=[src_ids_state, tar_ids_state],
                               outputs=[gr.Textbox(label="Landmark Setup Info")])

    run_button.click(fn=register_and_visualize_with_zip,
                     inputs=[target_input, source_input, output_dir, method_dropdown,
                             fricp_useinit, fricp_fileinit, fricp_maxiter, fricp_stop,
                             spare_iters, spare_stopcoarse, spare_stopfine],
                     outputs=[result_plot, download_button])

    rerun_button.click(fn=reuse_last_result_as_source,
                       inputs=[target_input, output_dir],
                       outputs=[source_input, upload_plot])

if __name__ == "__main__":
    demo.launch(server_name="0.0.0.0", server_port=7860)


