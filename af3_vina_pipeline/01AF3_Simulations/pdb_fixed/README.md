this contains the exact raw pdb files, but with the below required modifications.

If the protein says UNK:

these proteins were inferred locally on ICER, but with UNK residues in outputs
which errors the solvate script. As advised by Duncan, I just remove all UNK lines.


If the protein says CHAI:

these proteins where inferred locally on ICER, but gave "ERROR) Found 1 water atoms near the solute!" 
when trying to solvate them (reason unclear). Rerunning them on Chai-1 fixed the issue (Unsure
whether or run running on AF3 Server would have fixed.)

A list of these PDBs are saved in the 'errors' dir

As for specfically:
Human/A0A6B9RBU6 
DouglasFir/P85912
Both have X Amino Acid in chain, but AF3 local could not process them, so I resulted to CHAI, 
then cut out UNK manually


Species: Arabidopsis
  a0a178wen7.pdb: CHAI
  a0a654eb39.pdb: CHAI
  a0a654ebe9.pdb: CHAI
  q42266.pdb: UNK

Species: Human
  h0yd51.pdb: UNK
  f4mhi8.pdb: UNK
  a0a6b9rbu6.pdb: UNK
  a0a059qfd5.pdb: UNK
  h0yi80.pdb: UNK
  h7c4v4.pdb: UNK
  a0a059rhp1.pdb: UNK
  h0y7e5.pdb: UNK
  m0r351.pdb: UNK
  a0a343gqc8.pdb: UNK
  h3bqw1.pdb: UNK
  a0a3b3isn8.pdb: UNK
  h0yea5.pdb: UNK
  a0a286yf78.pdb: UNK
  e5dyd7.pdb: UNK
  a0a514yrt5.pdb: UNK
  a0a097q3h9.pdb: UNK
  h0yfr9.pdb: UNK
  a0a8f4gcf3.pdb: UNK
  a0a346m014.pdb: UNK
  h0ydg5.pdb: UNK
  a0a649yes9.pdb: UNK
  h0ycm4.pdb: UNK
  k7el57.pdb: UNK

Species: DouglasFir
  q14es8.pdb: UNK
  b8rin5.pdb: UNK
  p85912.pdb: UNK
  p85915.pdb: UNK
  p85909.pdb: UNK
  q85uy0.pdb: UNK
  c6752.pdb: CHAI

Species: Eucalyptus
  a0a059by44.pdb: CHAI
  a0a059axs5.pdb: CHAI
