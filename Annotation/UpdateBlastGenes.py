from multiprocessing import Pool, cpu_count
import time

# Load map file into a dictionary
map_file = "ficifolium_complete_deleted.map"
gff_file = "ficifolium.renamed.gff3"
output_file = "ficifolium.renamedBlast.gff3"

# Read map file and store replacements in a dictionary
id_map = {}
with open(map_file, 'r') as mf:
    for line in mf:
        old_id, new_id = line.strip().split()
        id_map[old_id] = new_id

# Function to replace IDs in a chunk of text without regex
def replace_ids_in_chunk(chunk_info):
    idx, chunk = chunk_info  # Unpack index and chunk text
    start_time = time.time()
    print(f"Processing chunk {idx}...")

    # Split by line to handle individual replacements within attributes
    updated_lines = []
    for line in chunk.splitlines():
        # Replace each old_id with new_id using .replace()
        for old_id, new_id in id_map.items():
            line = line.replace(old_id, new_id)
        updated_lines.append(line)
    
    # Join lines back after processing
    updated_chunk = "\n".join(updated_lines)
    print(f"Chunk {idx} processed in {time.time() - start_time:.2f} seconds.")
    return updated_chunk

# Read the entire GFF3 file into memory and split it into chunks for efficient parallel processing
with open(gff_file, 'r') as gf:
    gff_content = gf.read()

# Adjust chunk size based on file size and CPU count
num_chunks = cpu_count()  # Fewer chunks for efficient processing
chunk_size = len(gff_content) // num_chunks
chunks = [(i, gff_content[i:i + chunk_size]) for i in range(0, len(gff_content), chunk_size)]

# Process chunks in parallel
with Pool(cpu_count()) as pool:
    updated_chunks = pool.map(replace_ids_in_chunk, chunks)

# Write the updated content to the output file
with open(output_file, 'w') as out_f:
    out_f.write("".join(updated_chunks))

print(f"Updated GFF3 file saved to {output_file}")
