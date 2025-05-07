import pickle
import os

import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
from st_aggrid import AgGrid, GridOptionsBuilder, GridUpdateMode
import plotly.express as px


st.set_page_config(layout="wide", page_title="TF Atlas")
tab1, tab2 = st.tabs(["🧬 Benchmark Metadata", "HT-SELEX"])

# Load your DataFrame
@st.cache_data
def load_data():
    return pd.read_csv('/Users/andrewsg/Desktop/bioqc-results.txt', sep='\t', header=0)

with tab1:
    df = load_data()
    df['pearson'] = np.nan_to_num(df['pearson'])
    df['frac'] = np.nan_to_num(df['frac'])

    # Ensure 'notes' column exists
    if 'notes' not in df.columns:
        df['notes'] = ""
    
    
    # Fake sidebar layout
    sidebar_col, col1, col2 = st.columns([1, 4, 1])

    with sidebar_col:
        st.markdown("""
            <div style='background-color:#f9f9f9;padding:10px;border-radius:5px;'>
            <h4>🔍 Filter Options</h4>
        """, unsafe_allow_html=True)

        unique_tfs = df['tf'].dropna().unique()
        unique_cells = sorted(df['cell'].dropna().unique().tolist())
        unique_dbs = df['db'].dropna().unique()

        selected_dbs = st.multiselect("Select Database(s):", options=unique_dbs)
        selected_cells = st.multiselect("Select Cell(s):", options=unique_cells)
        selected_tfs = st.multiselect("Select Transcription Factor(s):", options=unique_tfs)

        st.markdown("""</div>""", unsafe_allow_html=True)

    # Filter the DataFrame
    filtered_df = df.copy()
    if selected_tfs:
        filtered_df = filtered_df[filtered_df['tf'].isin(selected_tfs)]
    if selected_cells:
        filtered_df = filtered_df[filtered_df['cell'].isin(selected_cells)]
    if selected_dbs:
        filtered_df = filtered_df[filtered_df['db'].isin(selected_dbs)]

    with sidebar_col:
        if selected_tfs:
            for tf in selected_tfs:
                motif_pkl = f'/Users/andrewsg/Desktop/canonical_motifs/{tf}.pkl'
                if os.path.exists(motif_pkl):
                    with open(motif_pkl, 'rb') as f:
                        fig = pickle.load(f)
                    st.pyplot(fig)

        st.markdown("---")
        st.markdown("### Thresholds")

        frac_range = st.slider(
            "Fraction of Peaks with Motif",
            min_value=0.0, max_value=1.0, value=(0.5, 1.0), step=0.01
        )

        pearson_range = st.slider(
            "Canonical Pearson Correlation",
            min_value=0.0, max_value=1.0, value=(0.9, 1.0), step=0.01
        )


        # Add QC pass column
        filtered_df['QC pass'] = (
            (filtered_df['pearson'] >= pearson_range[0]) &
            (filtered_df['pearson'] <= pearson_range[1]) &
            (filtered_df['frac'] >= frac_range[0]) &
            (filtered_df['frac'] <= frac_range[1])
        )

        # Apply pearson/frac threshold filters
        filtered_df = filtered_df[
            (filtered_df['pearson'] >= pearson_range[0]) &
            (filtered_df['pearson'] <= pearson_range[1]) &
            (filtered_df['frac'] >= frac_range[0]) &
            (filtered_df['frac'] <= frac_range[1])
        ]



    with col1:
        st.subheader("TF ChIP-seq Metadata")
        # pie_col1, pie_col2 = st.columns(2)

        # with pie_col1:
        #     db_counts = filtered_df['db'].value_counts()
        #     total = db_counts.sum()
        #     fig1, ax1 = plt.subplots()
        #     wedges, texts, autotexts = ax1.pie(
        #         db_counts,
        #         labels=db_counts.index,
        #         autopct=lambda pct: f"{int(pct/100*total)} ({pct:.1f}%)",
        #         startangle=90
        #     )
        #     ax1.set_title("Database Composition")
        #     st.pyplot(fig1)

        # if not selected_tfs:
        #     with pie_col2:
        #         unique_tf_count = filtered_df['tf'].nunique()
        #         remaining = 1639 - unique_tf_count
        #         fig2, ax2 = plt.subplots()
        #         ax2.pie([unique_tf_count, remaining],
        #                 labels=[f"Observed ({unique_tf_count})", f"Unobserved ({remaining})"],
        #                 autopct='%1.1f%%',
        #                 startangle=90)
        #         ax2.set_title("TF Coverage (out of 1,639)")
        #         st.pyplot(fig2)
        
        gb = GridOptionsBuilder.from_dataframe(filtered_df[['accession', 'db', 'tf', 'cell', 'pearson', 'frac', 'QC pass', 'notes']])
        gb.configure_selection(selection_mode="single", use_checkbox=True)
        gb.configure_column('notes', editable=True)
        grid_options = gb.build()

        response = AgGrid(
            filtered_df[['accession', 'db', 'tf', 'cell', 'pearson', 'frac', 'QC pass', 'notes']],
            gridOptions=grid_options,
            update_mode=GridUpdateMode.MODEL_CHANGED,
            theme="streamlit",
            enable_enterprise_modules=False,
            height=400,
            fit_columns_on_grid_load=True
        )

        # if response['data'] is not None:
        #     updated_data = pd.DataFrame(response['data'])

        #     # Reinsert any edited notes into the original df (identified by 'accession')
        #     for idx, row in updated_data.iterrows():
        #         accession = row['accession']
        #         new_note = row['notes']
        #         df.loc[df['accession'] == accession, 'notes'] = new_note

        if response["selected_rows"] is not None and len(response["selected_rows"]) > 0:
            st.write("You selected:", response.selected_data.accession[0])

        st.markdown("---")
        st.subheader("Motif Quality Metrics Distribution")
        
        # Count QC passing datasets
        qc_pass_df = filtered_df[filtered_df['QC pass']]
        num_passing = qc_pass_df.shape[0]
        total = df.shape[0]
        percent_passing = (num_passing / total) * 100 if total > 0 else 0
        unique_tfs_passing = qc_pass_df['tf'].nunique()

        # Display summary
        st.markdown(f"**✅ {num_passing:,} of {total:,} datasets pass QC ({percent_passing:.1f}%)**")
        st.markdown(f"**🧬 {unique_tfs_passing:,} unique TFs in passing datasets**")

        hist_col1, hist_col2, hist_col3 = st.columns(3)

        # Prepare masks
        frac_mask = (df['frac'] >= frac_range[0]) & (df['frac'] <= frac_range[1])
        pearson_mask = (df['pearson'] >= pearson_range[0]) & (df['pearson'] <= pearson_range[1])
        combined_mask = frac_mask & pearson_mask

        # Histogram: Fraction of Peaks
        with hist_col1:
            fig_frac, ax_frac = plt.subplots()
            all_vals = df['frac'].dropna()
            ax_frac.hist(all_vals, bins=30, color='skyblue', edgecolor='black')

            # Shade outside threshold
            ax_frac.axvspan(0, frac_range[0], color='gray', alpha=0.5)
            ax_frac.axvspan(frac_range[1], 1, color='gray', alpha=0.5)

            ax_frac.axvline(x=frac_range[0], linestyle='--', color='black')
            ax_frac.axvline(x=frac_range[1], linestyle='--', color='black')
            ax_frac.set_title("Fraction of Peaks with Motif")
            ax_frac.set_xlabel("Fraction")
            ax_frac.set_ylabel("Count")
            ax_frac.set_xlim([0, 1])
            st.pyplot(fig_frac)

        # Histogram: Pearson
        with hist_col2:
            fig_pearson, ax_pearson = plt.subplots()
            all_vals = df['pearson'].dropna()
            ax_pearson.hist(all_vals, bins=30, color='salmon', edgecolor='black')

            # Shade outside threshold
            ax_pearson.axvspan(0, pearson_range[0], color='gray', alpha=0.5)
            ax_pearson.axvspan(pearson_range[1], 1, color='gray', alpha=0.5)

            ax_pearson.axvline(x=pearson_range[0], linestyle='--', color='black')
            ax_pearson.axvline(x=pearson_range[1], linestyle='--', color='black')
            ax_pearson.set_title("Canonical Pearson Correlation")
            ax_pearson.set_xlabel("Pearson r")
            ax_pearson.set_ylabel("Count")
            ax_pearson.set_xlim([0, 1])
            st.pyplot(fig_pearson)

        # Scatter plot
        with hist_col3:
            fig_scatter, ax_scatter = plt.subplots()

            # Passing points
            ax_scatter.scatter(
                df.loc[combined_mask, 'frac'],
                df.loc[combined_mask, 'pearson'],
                color='purple', alpha=0.6, s=10, label='Pass'
            )

            # Failing points
            ax_scatter.scatter(
                df.loc[~combined_mask, 'frac'],
                df.loc[~combined_mask, 'pearson'],
                color='gray', alpha=0.3, s=10, label='Fail'
            )

            # Threshold lines
            ax_scatter.axvline(x=frac_range[0], linestyle='--', color='black')
            ax_scatter.axvline(x=frac_range[1], linestyle='--', color='black')
            ax_scatter.axhline(y=pearson_range[0], linestyle='--', color='black')
            ax_scatter.axhline(y=pearson_range[1], linestyle='--', color='black')

            ax_scatter.set_title("Frac vs Pearson")
            ax_scatter.set_xlabel("Fraction of Peaks with Motif")
            ax_scatter.set_ylabel("Canonical Pearson Correlation")
            ax_scatter.set_xlim([0, 1])
            ax_scatter.set_ylim([0, 1])
            ax_scatter.legend()
            st.pyplot(fig_scatter)





    with col2:
        st.subheader("Motifs")

        selected_motif = st.session_state.get("selected_motif", None)

        with st.container(height=1200):
            if response["selected_rows"] is not None and len(response["selected_rows"]) > 0:
                with open(f'/Users/andrewsg/Desktop/chip_pkl/{response.selected_data.accession[0]}.pkl', 'rb') as f:
                    motifs_dict = pickle.load(f)

                sort_by = st.radio(
                    "Sort motifs by:",
                    options=["Canonical Similarity", "Fraction of Peaks"],
                    horizontal=True
                )

                # Sort motifs_dict keys based on the selected sorting method
                if sort_by == "Canonical Similarity":
                    sorted_items = sorted(motifs_dict.items(), key=lambda x: x[1].get("canonical_pearson", 0), reverse=True)
                elif sort_by == "Fraction of Peaks":
                    sorted_items = sorted(motifs_dict.items(), key=lambda x: x[1].get("frac", 0), reverse=True)
                else:
                    sorted_items = motifs_dict.items()

                for i, (motif_idx, motif_data) in enumerate(sorted_items):
                    fig, ax = plt.subplots(figsize=(1, 0.25), dpi=1200, tight_layout=True)
                    n = motif_data['n']
                    auc = motif_data['auc']
                    frac = motif_data['frac']
                    ppm = pd.DataFrame(motif_data['ppm'], columns=['A', 'C', 'G', 'T'])

                    logomaker.Logo(ppm * np.log2((ppm + .01) / 0.25), ax=ax)

                    ax.set_ylim([0, 2])
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                    ax.spines['bottom'].set_visible(False)
                    ax.spines['left'].set_linewidth(0.25)
                    ax.set_yticks([0, 1, 2])
                    ax.tick_params(axis='y', direction='in', labelsize=4, width=0.25)
                    ax.set_xticks([])
                    #ax.set_title(f'n={n:,} ({100*frac:.2f}%); auROC={auc:.2f}\n canonical pearson = {motif_data['canonical_pearson']:.2f}', fontsize=4)

                    plt.tight_layout(pad=0.0)
                    st.pyplot(fig)
                    st.write(f'n={n:,} ({100*frac:.2f}%); auROC={auc:.2f};\n canonical pearson = {motif_data['canonical_pearson']:.2f}')

                    col_c, col_s, col_d, col_a = st.columns(4)

                    centrality_state_key = f"centrality_state_{motif_idx}"  # for session state
                    centrality_button_key = f"centrality_btn_{motif_idx}"   # for Streamlit button widget

                    with col_c:
                        if st.button("Centrality", key=centrality_button_key):
                            current = st.session_state.get(centrality_state_key, False)
                            st.session_state[centrality_state_key] = not current


                    with col_s:
                        st.markdown(
                            f"<button style='font-size:10px;padding:2px 6px;border-radius:4px;'>Conservation</button>",
                            unsafe_allow_html=True
                        )

                    with col_d:
                        st.markdown(
                            f"<button style='font-size:10px;padding:2px 6px;border-radius:4px;'>Download</button>",
                            unsafe_allow_html=True
                        )

                    with col_a:
                        st.markdown(
                            f"<button style='font-size:10px;padding:2px 6px;border-radius:4px;'>Add to Cart</button>",
                            unsafe_allow_html=True
                        )
                    
                    if st.session_state.get(centrality_state_key, False):
                        st.markdown("**Centrality Distribution (placeholder)**")
                        fig_centrality, ax_centrality = plt.subplots()
                        ax_centrality.hist(np.random.normal(loc=0, scale=1, size=100), bins=20, color='lightblue', edgecolor='black')
                        ax_centrality.set_title("Centrality Placeholder Histogram")
                        st.pyplot(fig_centrality)



