import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import plotly.express as px
import textwrap
from concurrent.futures import ProcessPoolExecutor, as_completed

def load_goea_results(path, dictionary):
    for file in path.glob("*_GOEA_results.txt"):
            goea_results_df = pd.read_csv(file, sep='\t')
            file_name = file.stem.replace("_GOEA_results", "").replace("_", " ")
            dictionary[file_name] = goea_results_df

def wrap_text(text, width=30):
    """Wraps text to a specified width for better readability."""
    return "<br>".join(textwrap.wrap(text, width))

def create_enrichment_barplots(goea_results_df, taxonomic_level, name, outdir):

    min_value = np.nextafter(0, 1)
    goea_results_df['-log10(Stat. Sig.)'] = -np.log10(
        goea_results_df['Statistical significance'].replace(0, min_value)
    )
    goea_results_df['Formatted Statistical significance'] = goea_results_df['Statistical significance'].apply(lambda x: f"{x:.2e}")
    goea_results_df['Odds Ratio'] = np.exp(goea_results_df['log(Odds Ratio)'])
    goea_results_df['Formatted Odds Ratio'] = goea_results_df['Odds Ratio'].apply(lambda x: f"{x:.2f}")

    for namespace in ["Biological Process", "Molecular Function", "Cellular Component"]:
        
        namespace_df = goea_results_df[(goea_results_df["Namespace"] == namespace)]
        if namespace_df.empty:
            print(f"No valid {namespace} results for {taxonomic_level}.")
            continue  

        sorted_namespace_df = namespace_df.sort_values(by="Statistical significance", ascending=True).head(20)
        sorted_namespace_df = sorted_namespace_df.sort_values(by="log(Odds Ratio)", ascending=True)
        wrapped_names = sorted_namespace_df['Name'].apply(lambda x: wrap_text(x, width=30))
        
        if taxonomic_level == "taxonomy":
            sorted_namespace_df['Coverage'] = (sorted_namespace_df['N with GO term'].astype(str) + "/" +
                                               sorted_namespace_df['N in taxonomy'].astype(str))
            sorted_namespace_df['Formatted Percentage'] = sorted_namespace_df['Species Percentage'].apply(lambda x: f"{x:.1f}%")

            custom_cols = ['Name',
                           'GO Term ID',
                           'Formatted Statistical significance',
                           '-log10(Stat. Sig.)',
                           'Formatted Percentage',
                           'Coverage',
                           'Heterogeneity',
                           'Formatted Odds Ratio']
            hover_template = (
                "<b>%{customdata[0]}</b><br>" +
                "<b>GO ID:</b> %{customdata[1]}<br>" +
                "<b>log(Odds Ratio):</b> %{x:.2f}<br>" +
                "<b>Odds Ratio:</b> %{customdata[7]}<br>" +
                "<b>Stat. Sig.:</b> %{customdata[2]}<br>" +
                "<b>-log10(Stat. Sig.):</b> %{customdata[3]:.2f}<br>" +
                "<b>Taxonomic coverage:</b> %{customdata[5]} (%{customdata[4]})<br>" +
                "<b>Heterogeneity:</b> %{customdata[6]:.2f}"
            )
        elif taxonomic_level == "species":
            custom_cols = ['Name',
                           'GO Term ID',
                           'Formatted Statistical significance',
                           '-log10(Stat. Sig.)',
                           'N Study with term',
                           'N Study without term',
                           'N Background with term',
                           'N Background without term',
                           'Formatted Odds Ratio']
            hover_template = (
                "<b>%{customdata[0]}</b><br>" +
                "<b>GO ID:</b> %{customdata[1]}<br>" +
                "<b>log(Odds Ratio):</b> %{x:.2f}<br>" +
                "<b>Odds Ratio:</b> %{customdata[8]}<br>" +
                "<b>Stat. Sig.:</b> %{customdata[2]}<br>" +
                "<b>-log10(Stat. Sig.):</b> %{customdata[3]:.2f}<br>" +
                "<b>Number of Study proteins with term:</b> %{customdata[4]}<br>" +
                "<b>Number of Study proteins without term:</b> %{customdata[5]}<br>" +
                "<b>Number of Background proteins with term:</b> %{customdata[6]}<br>" +
                "<b>Number of Background protein without term:</b> %{customdata[7]}"
            )
        else:
            raise ValueError("Invalid taxonomic_level. Choose 'taxonomy' or 'species'.")

        fig = px.bar(
            sorted_namespace_df,
            y=wrapped_names,
            x='log(Odds Ratio)',
            color='-log10(Stat. Sig.)',
            color_continuous_scale="Cividis",
            labels={'log(Odds Ratio)': 'log(Odds Ratio) (Higher is better)', 'y': ''},
            custom_data = custom_cols,
            width=600,
            height=900,
        )
        
        fig.update_traces(
            hovertemplate=hover_template
        )
        
        fig.update_coloraxes(
            colorbar_title="-log10(Stat. Sig.)",
            colorbar_tickfont=dict(size=10),
            colorbar_title_font=dict(size=12),
            colorbar_len=0.2,
            colorbar_thickness=15,
            colorbar_x=1,
            colorbar_y=0.91,
            colorbar_yanchor="middle"
        )
        
        fig.update_layout(
            margin=dict(l=10, r=10, t=30, b=10),
            showlegend=False,
            dragmode=False,
            xaxis=dict(
                title_font=dict(size=12),
                tickfont=dict(size=10)
            ),
            yaxis=dict(
                tickfont=dict(size=12)
            ),
            bargap=0.4,
            yaxis_showgrid=True
        )

        fig.write_html(f"{outdir}/{namespace}/{name}.html")
        fig.write_image(f"{outdir}/{namespace}/{name}.png")
def process_species(task):
    """Wrapper to call the plotting function for one species."""
    name, goea_results_df, taxonomic_level, outdir = task
    create_enrichment_barplots(goea_results_df, taxonomic_level, name, outdir)

def main():
    parser = argparse.ArgumentParser(description="Example script")
    parser.add_argument("-i", help="Directory containing the GOEA results.", required=True)
    args = parser.parse_args()

    goea_results_dir = Path(args.i)
    goea_results_dict = dict()

    if (goea_results_dir / "single_taxon_results").exists():
        results_dir = goea_results_dir / "single_taxon_results"
        taxonomic_level = "species"
    else:
        print("No valid result directories found.")
        return

    load_goea_results(results_dir, goea_results_dict)
    barplots_dir = results_dir / "plots" / "barplots"
    barplots_dir.mkdir(parents=True, exist_ok=True)

    for namespace in ["Biological Process", "Molecular Function", "Cellular Component"]:
        (barplots_dir / namespace).mkdir(parents=True, exist_ok=True)

    tasks = [
        (name, goea_results_df, taxonomic_level, barplots_dir)
        for name, goea_results_df in goea_results_dict.items()
    ]

    with ProcessPoolExecutor() as executor:
        futures = [executor.submit(process_species, task) for task in tasks]
        for future in as_completed(futures):
            try:
                future.result()
            except Exception as e:
                print(f"Error processing species: {e}")

if __name__ == "__main__":
    main()