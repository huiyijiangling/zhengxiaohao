# Preprocessing of manufacturers' annotation files to make them smaller and more manageable

options(stringsAsFactors = FALSE);

tab = read.delim("../Raw-Large/Agilent014850_D_AA_20070207-Tomida.txt");

tab2 = tab[, c(1:10)];

write.table(tab2, file = bzfile("../Shortened/Agilent014850_D_AA_20070207-Tomida-shortened.txt.bz2"),
            sep = "\t", quote = FALSE, row.names = FALSE);

apply(tab2, 2, mode)

