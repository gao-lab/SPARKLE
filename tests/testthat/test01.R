# load(file = "data/cwas.test.data.rda")
# load(file = "data/cwas.test.data.single.model.rda")

#cwas_gene_mediate_cell_analysis_singlecore(cwas.data = cwas.test.data,Best_model = cwas.test.data.single.model.best,X.cell = "MTF",Mediation = "TGFB1")

#mytest.result1 <- readRDS("mytest.result1/CWAS_resultmytest.result12024-05-22105225.rds")

# test_that("save correctly", {
#
#   result1 <- cwas_allmodel_cal(cwas.test.data)
#   result1$`MYO-F`$Models_info$id
#   expect_equal(result1$`MYO-F`$Models_info$id, cwas.test.data.single.model$`MYO-F`$Models_info$id )
#
# })

# test_that("mediation calculation", {
#
#   result1 <- mediation_analysis_auto(cwas.data = cwas.test.data,Best_model = cwas.test.data.single.model.best,
#                                      X.cell = "MTF",X.gene = "TGFB1",Moderation.cell ="APF1",Moderation.gene = "TGFB2",nsim = 100,seed = 1,CI = "boot" )
#
#   #sf.rds(result1,"mytest.result1")
#  # expect_equal(result1$results[[1]]$mediation$pval, mytest.result1$results[[1]]$mediation$pval )
#
# })

# test_that("mediation calculation", {
#
#
#   result1 <- cwas_mediation_analysis(cwas.data = cwas.test.data,Best_model = cwas.test.data.single.model.best,
#                                      X.cell = "MTF",X.gene = "TGFB1",mediate.cell = "APF1",
#                                      mediate.gene =  "TGFB2",method = "gene_mediate_cell",fix_effect = "Group",random_effect = "Subgroup",num_cores = 2 )
#
#   #sf.rds(result1,"mytest.result1")
#   # expect_equal(result1$results[[1]]$mediation$pval, mytest.result1$results[[1]]$mediation$pval )
#
# })

test_that("mediation calculation", {


  result1 <- cwas_moderation_analysis(cwas.data = cwas.test.data,Best_model = cwas.test.data.single.model.best,
                                     X.cell = "MTF",X.gene = "TGFB1",Moderation.cell = "APF1",
                                     Moderation.gene =  "TGFB2",method = "cell_moderation_gene",num_cores = 2 )

  #sf.rds(result1,"mytest.result1")
  # expect_equal(result1$results[[1]]$mediation$pval, mytest.result1$results[[1]]$mediation$pval )

})


###################################
