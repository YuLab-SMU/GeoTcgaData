test_that("can parse example id_conversion_vector", {
    # skip_on_bioc()
    result <- id_conversion_vector(
        "symbol", "ensembl_gene_id",
        c("A2ML1", "A2ML1-AS1", "A4GALT", "A12M1", "AAAS"))
    expect_equal(colnames(result) , c("from", "to"))
})