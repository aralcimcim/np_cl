# chunkLearner_IBP
chunk learner with nonparametric Bayesian modeling using the Indian Buffet Process

- codes for running the model:
    wood_batch.m                    script for parameter and input definition and running the model that calculates the posteriors
    wood_ibp_learning_frontend.m    function to call the Gibbs sampler for Y and Z
    wood_make_gibbs_y_spatial.m     Gibbs sampler for Y
    wood_make_gibbs_z_spatial.m     Gibbs sampler for Z
    calc_px_training.m              calculation of likelihood P(X_t|Z, Y_:t)
    calc_pv_training.m              calculation of likelihood P(V|Z, Y)
    calc_pv_trial.m                 trial-wise auxiliary calculaion for calc_pv_training

- codes that might be useful later for scene probabilities and choice probabilities
    humanScript.m                   script for running the model for whole experiments to return choice probabilities
    humantest.m                     function for calculating the choice probabilities from the posteriors

    pScene_spatial.m                probability of a scene, given the posterior Z and R 
                                    P(XScene) =
                                        = sum_{all output datum: for all steps} P(XScene | Z, R) * P(Z, R) =
                                        = sum pScene_ZR(XScene, Z, R..) * 1/stepNo
        pScene_ZR_spatial.m         auxiliary fn for pScene_spatial

    pScene_spatial_Cpost.m          probability of a scene, given the posterior Z and R using the C posterior 
        pScene_ZR_spatial_Cpost.m   auxiliary fn for pScene_spatial_Cpost
        calc_pv_training_Cpost.m    auxiliary fn for pScene_ZR_spatial_Cpost
        calc_pv_trial_Cpost.m       auxiliary fn for calc_pv_training_Cpost


- inputs from the 2001 experiment and a truncated version are in the inputs folder
- model doc is IBP_modelDoc.pdf
- original paper for a parametric version of this model is Orban_PNAS_2008
- the IBP is derived from Wood_IBP_2012.pdf, with code ibp_wood.zip

- the python folder contains the python implementation from Aral and its commented alternative (comments from Tünde)
- the nonSpatial folder contains a model setup where the spatial information is not exploited hence the model runs much faster 

