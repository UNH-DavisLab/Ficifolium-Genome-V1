def extract_gene_locations(gff_file_path):
    gene_locations = set()  # To store unique gene entries

    with open(gff_file_path, 'r') as file:
        for line in file:
            # Only process lines with the 'gene' feature type
            if '\tgene\t' in line:
                columns = line.strip().split('\t')
                chromosome = columns[0]
                start = columns[3]
                stop = columns[4]

                # Extract gene ID from the attributes column
                attributes = columns[8]
                gene_id = None
                for attribute in attributes.split(';'):
                    if attribute.startswith('ID='):
                        gene_id = attribute.split('=')[1]
                        break

                # Add the unique entry to the set if gene_id was found
                if gene_id:
                    gene_locations.add((chromosome, gene_id, start, stop))

    # Print results in the specified format
    for chromosome, gene_id, start, stop in sorted(gene_locations):
        print(f"{chromosome}\t{gene_id}\t{start}\t{stop}")

gff_file_path = 'ficifolium_final_annotation_sorted.gff3'
extract_gene_locations(gff_file_path)
