from Bio import Entrez, SeqIO


def get_numbers_list(file):
    file = open(file, "r")
    mydata = file.read()
    dataintolist = mydata.split("\n")
    print(dataintolist)
    file.close()
    return dataintolist


def get_sequences(numbers):
    Entrez.email = "beatrizf.2000@gmail.com"
    for n in numbers:
        handle = Entrez.efetch(db="nucleotides", rettype="gb", retmode="text", id=n)
        seq_record = SeqIO.read(handle, "db")
        print(seq_record)
        SeqIO.write(seq_record,'C:\Users\beatr\Desktop\Universidade\projeto\genoma' + seq_record.id[:-2] + ".gb", "gb")


def gbk_to_fasta(gbk_filename, new_fasta_name):
    input_handle = open(gbk_filename, "r")
    output_handle = open(new_fasta_name, "w")

    for seq_record in SeqIO.parse(input_handle, "genbank"):
        print("Dealing with GenBank record %s" % seq_record.id)
        output_handle.write(">%s %s\n%s\n" % (seq_record.id,seq_record.description,seq_record.seq))

    output_handle.close()
    input_handle.close()

gbk_to_fasta("C:\Users\beatr\Desktop\Universidade\projeto\genoma\NC_001341.gb", "C:\Users\beatr\Desktop\Universidade\projeto\genoma\NC_001341.fasta")