# ccbr1046 Project 
CRISPRESSO Annotation

1_Align_Code.R Parses through the CRSIPResso output to run a local alignment and prepares a file in AVINPUT format for ANNOVAR.
2_runanno.sh Uses the Output from the previous step to call ANNOVAR with specified annotations (CLINVAR, SIFT, POLYPHEN, etc)
3_ANNOVAR_Parse.R Parses the ANNOVAR Output files (multianno.txt)
