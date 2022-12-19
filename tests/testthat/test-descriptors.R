
test_that('smartsSearch errors if > 1 SMILES specified',{
  expect_error(smartsSearch(amino_acids$SMILES,"[OX2H]"))
})