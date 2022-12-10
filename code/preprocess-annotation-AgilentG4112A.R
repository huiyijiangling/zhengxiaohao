# Preprocess the file annotation-AgilentG4112A-FromEBI.txt.bz2 to remove extra columns and rename the
# remaining ones in a more descriptive form

library(WGCNA)
options(stringsAsFactors = FALSE);
annot0 = read.delim(bzfile("../annotation-AgilentG4112A-FromEBI.txt.bz2"));

drop = c("Reporter.Name", "Comment.AEReporterName.", "Reporter.Database.Entry.agilent_control.EQC1.0.", 
         "Reporter.Sequence", "Reporter.Comment", "Composite.Element.Name",
         "Composite.Element.Database.Entry.agilent_control.EQC1.0.", "Composite.Element.Comment");

dropInd = match(drop, colnames(annot0));

annot1 = annot0[, -dropInd];

colnames(annot1) = gsub("Composite.Element.Database.Entry.", "", colnames(annot1), fixed = TRUE);
colnames(annot1) = gsub("\\.$", "", colnames(annot1), fixed = FALSE);

write.table(annot1, file = bzfile("../annotation-AgilentG4112A-FromEBI-processed.txt.bz2"),
            quote = FALSE, sep = "\t", row.names = FALSE);
