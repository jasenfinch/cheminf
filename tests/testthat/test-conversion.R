
test_that('conversion from SMILES to INCHI works',{
  skip_if(Sys.info()["sysname"] == 'Windows','Skipping test on Windows.')
  inchi <- convert(amino_acids$SMILES[1],'SMILES','INCHI')
  expect_true(identical(inchi,"InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H2,1H3,(H,5,6)/t2-/m0/s1"))
})

test_that('smilesToMF errors if > 1 SMILES specified',{
  expect_error(smilesToMF(amino_acids$SMILES))
})


test_that('smilesToAccurateMass errors if > 1 SMILES specified',{
  expect_error(smilesToAccurateMass(amino_acids$SMILES))
})
