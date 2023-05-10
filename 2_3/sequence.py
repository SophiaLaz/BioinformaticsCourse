def generate_randon_seq(n: int) -> str:
    import random
    return ''.join(random.choice('ATGC') for _ in range(n))


def k_mers(seq: str, length: int) -> list:
    return [seq[i:i + length] for i in range(len(seq) - length + 1)]


def generate_reads(seq: str, length: int, count: int = 0) -> list:
    import random
    len_str = len(seq)
    count = len_str - length + 1 if count == 0 else count
    indexes = [random.randint(0, len_str - length) for _ in range(count)]
    return [seq[ind:ind + length] for ind in indexes]


def read_fasta(file_name: str) -> list:
    from Bio import SeqIO
    file_name = file_name if file_name[-6:] == '.fasta' else file_name + '.fasta'
    data = list(SeqIO.parse(file_name, 'fasta'))
    return [str(data[i].seq) for i in range(len(data))]


def create_fasta(seq: list, file_name: str):
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    file_name = file_name if file_name[-6:] == '.fasta' else file_name + '.fasta'
    records = [SeqRecord(Seq(seq[i]), id='', description='') for i in range(len(seq))]
    SeqIO.write(records, file_name, 'fasta')


def fastq_to_fasta(file_name: str) -> list:
    from Bio import SeqIO
    records = SeqIO.to_dict(SeqIO.parse(file_name, "fastq"))
    data_fasta = [records[records[i]].format("fasta") for i in range(len(records))]
    return [str(data_fasta[i].seq) for i in range(len(data_fasta))]


def read_fastq(file_name: str) -> list:
    from Bio import SeqIO
    file_name = file_name if file_name[-6:] == '.fastq' else file_name + '.fastq'
    data = list(SeqIO.parse(file_name, 'fastq'))
    return [str(data[i].seq) for i in range(len(data))]
