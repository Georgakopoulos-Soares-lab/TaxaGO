import pandas as pd 
from pathlib import Path
import urllib.request
import urllib.error
import urllib.parse

ftp_urls = [
    "ftp://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/proteome2taxid",
    "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes/README"
]

print("Starting download process...")

for ftp_url in ftp_urls:
    print("-" * 40)

    parsed_url = urllib.parse.urlparse(ftp_url)
    url_path_str = parsed_url.path

    local_filename_str = Path(url_path_str).name

    local_path = Path(local_filename_str)

    full_local_path = local_path.resolve()
    print(f"Processing URL: {ftp_url}")
    print(f"Target local file: {full_local_path}")
    print(f"Checking existence...")

    if not local_path.exists():
        print(f"-> File '{local_path.name}' not found locally. Attempting download...")
        urllib.request.urlretrieve(ftp_url, local_path)
        print(f"-> Download complete! File saved as '{local_path}'.")
        
    else:
        print(f"-> File '{local_path.name}' already exists. Skipping download.")

print("-" * 40)
print("Download process finished.")

header_line = "Proteome_ID\tTax_ID\tOSCODE\tSUPERREGNUM\t#(1)\t#(2)\t#(3)\tSpecies Name"
found_header = False
lines_written = 0
input_filename = "README"
output_filename = "processed_README"

with open(input_filename, 'r', encoding='utf-8') as infile, \
        open(output_filename, 'w', encoding='utf-8') as outfile:

    for line in infile:
        stripped_line = line.strip()

        if found_header:
            if stripped_line == "":
                break
            else:
                outfile.write(line)
                lines_written += 1
        elif stripped_line == header_line:
            found_header = True
            outfile.write(line)
            lines_written += 1

if not found_header:
    print(f"Error: Header line not found in '{input_filename}'")
elif lines_written == 1 and found_header:
        print(f"Header found, but no data lines followed before the end or a blank line in '{input_filename}'. Output file '{output_filename}' contains only the header.")
elif lines_written > 1:
    print(f"Successfully extracted data to '{output_filename}'")
