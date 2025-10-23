import os
import re

# Define base paths
source_root = "../SmallMolecules/cgenff_results"
destination_root = "../SmallMolecules/topologies"

# Make sure destination directory exists
os.makedirs(destination_root, exist_ok=True)

# Walk through all subdirectories and files in source_root
for dirpath, _, filenames in os.walk(source_root):
    for filename in filenames:
        if filename.endswith(".str"):
            file_path = os.path.join(dirpath, filename)

            # Read content of the .str file
            with open(file_path, 'r') as f:
                content = f.read()

            # Replace the word after RESI
            base_name = os.path.splitext(filename)[0]
            replacement = f"L{base_name}"
            content_modified = re.sub(r'(RESI\s+)(\S+)', r'\1' + replacement, content)

            # Save to topologies directory using same filename
            new_file_path = os.path.join(destination_root, filename)
            with open(new_file_path, 'w') as f:
                f.write(content_modified)

print("All .str files processed and saved to '../SmallMolecules/topologies'.")
