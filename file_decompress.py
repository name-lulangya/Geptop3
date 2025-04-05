import zipfile
import os
from pathlib import Path

def decompress_zips(zip_dir, target_dir):
    """
    Decompress all ZIP files in a directory and rename extracted files with .pkl extension
    
    Args:
        zip_dir (str): Directory path containing ZIP files
        target_dir (str): Target directory for decompressed files
    """
    # Create target directory if not exists
    Path(target_dir).mkdir(parents=True, exist_ok=True)
    
    # Process all ZIP files in the directory
    for zip_file in Path(zip_dir).glob("*.zip"):
        try:
            with zipfile.ZipFile(zip_file, 'r') as zipf:
                # Get first file in ZIP (assuming single file per archive)
                file_list = zipf.namelist()
                if not file_list:
                    print(f"Warning: {zip_file.name} is empty, skipped")
                    continue
                
                original_file = file_list[0]
                # Generate new filename with .pkl extension
                new_filename = f"{Path(original_file).stem}.pkl"
                target_path = Path(target_dir) / new_filename
                
                # Prevent overwriting existing files
                if target_path.exists():
                    print(f"Warning: {target_path} already exists, skipping")
                    continue
                
                # Extract and rename file
                with zipf.open(original_file) as src_file:
                    content = src_file.read()
                    with open(target_path, 'wb') as dst_file:
                        dst_file.write(content)
                print(f"Decompressed: {zip_file.name} → {new_filename}")
                
        except Exception as e:
            print(f"Error: Failed to decompress {zip_file.name} - {str(e)}")

# Usage example
root_dir = ""
decompress_zips(
    zip_dir = f"{root_dir}/CVFile_compress",  # Source directory with ZIP files
    target_dir = f"{root_dir}/CVFile"  # Target output directory
)