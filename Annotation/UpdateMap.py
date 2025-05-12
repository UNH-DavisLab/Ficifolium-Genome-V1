import re

gff_file = "ficifolium.gff3"
map_file = "ficifolium.map"
updated_map_file = "ficifolium_complete.map"

# Load the existing map file into a dictionary
mapped_ids = {}
with open(map_file, 'r') as mf:
    for line in mf:
        old_id, new_id = line.strip().split()
        mapped_ids[old_id] = new_id

# Extract all gene IDs from the GFF3 file
gff_ids = set()
with open(gff_file, 'r') as gf:
    for line in gf:
        # Match IDs starting with 'g' followed by 1-5 digits, .t, and a final digit
        matches = re.findall(r'\b(g\d{1,5}\.t\d)\b', line)
        gff_ids.update(matches)

# Find missing IDs
missing_ids = gff_ids - mapped_ids.keys()

# Generate new IDs for missing entries
new_entries = []

# Extract numeric portions from existing Cfic_ IDs, ignoring any suffixes
next_id_number = (
    max(
        int(re.search(r'Cfic_(\d+)', new_id).group(1))
        for new_id in mapped_ids.values()
        if re.search(r'Cfic_(\d+)', new_id)
    ) + 1
)

for missing_id in sorted(missing_ids):
    new_id = f"Cfic_{str(next_id_number).zfill(8)}"
    new_entries.append(f"{missing_id}\t{new_id}")
    next_id_number += 1

# Write the complete mapping file
with open(updated_map_file, 'w') as out_mf:
    # Copy existing mappings
    with open(map_file, 'r') as mf:
        out_mf.write(mf.read())
    # Add new entries
    for entry in new_entries:
        out_mf.write(entry + "\n")

print(f"Updated mapping file with missing IDs saved to {updated_map_file}")
