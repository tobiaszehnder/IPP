"""Create dummy data for testing."""

import os
import struct


def create_dummy_pwaln_bin(file_path: str):
    """
    Creates a minimal valid dummy.pwaln.bin file matching the binary format
    expected by Ipp::loadPwalns().

    The structure is:
    - version [uint8]
    - endianness_magic [uint16] (0xAFFE)
    - num_sp1 [uint8]
    - For each sp1:
        - sp1_name [null-terminated string]
        - sp1_genome_size [uint64]
        - num_sp2 [uint8]
        - For each sp2:
            - sp2_name [null-terminated string]
            - num_ref_chrom_entries [uint32]
            - For each ref chrom:
                - ref_chrom [uint32]
                - num_pwaln_entries [uint32]
                - For each pwaln entry:
                    - ref_start [uint32]
                    - qry_start [uint32]
                    - qry_chrom [uint32]
                    - length_and_strand [uint16]
    - num_chromosomes [uint32]
    - For each chromosome:
        - chrom_name [null-terminated string]

    This dummy will contain:
    - 1 species (sp1): "mm39" genome size 2000000
    - 1 query species (sp2): "hg38"
    - 1 reference chromosome (id=0)
    - 2 pwaln entries
    - 2 chromosomes: "chr1", "chr2"
    """
    version = 4
    endianness_magic = 0xAFFE
    sp1_name = b"mm39\x00"
    sp1_genome_size = 2000000
    sp2_name = b"hg38\x00"
    num_ref_chrom_entries = 1
    ref_chrom_id = 0
    num_pwaln_entries = 2

    # Define two dummy pwaln entries
    # Each entry is 16 bytes padded struct, but only 14 bytes packed in file
    # We'll pack fields separately (14 bytes):
    # ref_start (4), qry_start (4), qry_chrom (4), length_and_strand (2)
    pwaln_entries = [
        (100, 150, 0, 50),        # length_and_strand MSB=0 (forward strand)
        (300, 350, 1, 0x800A)     # length 10 with MSB=1 (reverse strand)
    ]

    chromosomes = [b"chr1\x00", b"chr2\x00"]

    with open(file_path, "wb") as f:
        # version
        f.write(struct.pack("B", version))
        # endianness magic (little endian)
        f.write(struct.pack("<H", endianness_magic))
        # num_sp1
        f.write(struct.pack("B", 1))

        # sp1_name
        f.write(sp1_name)
        # sp1_genome_size
        f.write(struct.pack("<Q", sp1_genome_size))
        # num_sp2
        f.write(struct.pack("B", 1))

        # sp2_name
        f.write(sp2_name)
        # num_ref_chrom_entries
        f.write(struct.pack("<I", num_ref_chrom_entries))

        # For each ref chrom
        # ref_chrom
        f.write(struct.pack("<I", ref_chrom_id))
        # num_pwaln_entries
        f.write(struct.pack("<I", num_pwaln_entries))

        # pwaln entries (packed, 14 bytes each)
        for ref_start, qry_start, qry_chrom, length_and_strand in pwaln_entries:
            # Pack as <IIIH (little-endian uint32, uint32, uint32, uint16)
            f.write(struct.pack(
                "<IIIH", ref_start, qry_start, qry_chrom,length_and_strand
            ))

        # num_chromosomes
        f.write(struct.pack("<I", len(chromosomes)))
        # chromosomes (null terminated strings)
        for chrom in chromosomes:
            f.write(chrom)

    print(f"Dummy pwaln.bin file created at {file_path}")

def create_dummy_bed_file(file_path: str):
    """
    Creates a dummy BED file with a few regions matching chromosomes
    present in the dummy pwaln.bin file.
    """
    bed_content = """chr1\t90\t160\tregion1
chr1\t280\t360\tregion2
chr2\t100\t200\tregion3
chr2\t400\t420\tregion4
"""
    with open(file_path, "w") as f:
        f.write(bed_content)
    print(f"Dummy BED file created at {file_path}")

def main():
    data_dir = "data"
    os.makedirs(data_dir, exist_ok=True)

    pwaln_path = os.path.join(data_dir, "dummy.pwaln.bin")
    bed_path = os.path.join(data_dir, "dummy_regions.bed")

    create_dummy_pwaln_bin(pwaln_path)
    create_dummy_bed_file(bed_path)

if __name__ == "__main__":
    main()