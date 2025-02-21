import plotly.express as px
from pathlib import Path
import argparse
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Example script")
    parser.add_argument("--input-dir", help="Directory containing the GOEA results.", required=True)
    # parser.add_argument("--output-dir", help="Directory to save the created figures", required=True)
    
    args = parser.parse_args()

    # Sample data
    data = {
        "Fruits": ["Apples", "Oranges", "Bananas", "Grapes"],
        "Amount": [10, 15, 7, 12]
    }

    # Create bar chart
    fig = px.bar(data, x="Fruits", y="Amount", title="Fruit Amounts")

    # # Save chart as a static image (e.g., PNG)
    # output_path = Path(args.input_dir) / "fruit_amounts.png"
    # fig.write_image(str(output_path))

    fig.show(renderer="svg")
    # print(f"Figure saved at: {output_path}")

if __name__ == "__main__":
    main()
