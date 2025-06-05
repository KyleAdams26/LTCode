This repository contains the code used in "Understanding Immune Dynamics in Liver Transplant Through Mathematical Modeling" by Bruner et. al. The following are instructions on how to use the MATLAB scripts to run your own sensitivity analysis and model simulations using your own model.

The following document is a walkthrough of how to use the MATLAB code to conduct a sensitivity analysis on your own model using the code from the Liver Transplant Rejection project.

1. Download all scripts
calculateQOI.m
calculateSobolIndices.m
compareQOIDistributions.m
generateRandomSubset.m
generateSobolSamples.m
makeSensitivityTables.m
odefun.m
plotModelSimulation.m
plotsVaryingInfluentialParams.m
setInitialConditions.m
setParameters.m
sobolMain.m
transformSobolSampleDistributions.m

2. Adjust odefun.m to your model.
In odefun.m, start with "% ###Step 1". The vector “c” represents your state variables, and you can adjust the length of this vector here if needed. If you model as x_1, … , x_n state variables, you will have n lines here of x_i = c(i) where you can name x_i something meaningful.
Next, see "% $$$Step 2". Here we have added the “pathways” of our model. Our pathways here are subsections of our equations. These pathways are made up of the state variable names you made in 2a and the parameters from the required structure “p” from the inputs of odefun.m. This “p” structure will have the names of our parameters along with their nominal values, and will be created in a separate file. Each parameter here should follow the example shown, in the form of p.parameterName.
Next, see "% ###Step 3". These should represent the dynamics/equations of your model. It is important these are in the same order as the state variables declared in step 2a. If in step 2a you declared x1 = c(1), then in this step, dy(1) must represent the derivative of x1 with respect to your independent variable. The number of these equations should also match the number of your state variables in step 2a. These dynamics should only refer to the parameters you have named in step 2c, and the variables you have named in step 2b. 
Lastly, adjust dcdt to contain as many equations as you made in step 2d.

3. Adjust setParameters.m
Here you should create a line for each parameter, where the name of the parameter matches the names you used in step 2c. For each line, set p.parameterName = parameterValue, as shown in the example.
This script will be used anytime you need to load parameters into another script.

4. Adjust setInitialConditions.m
Here is a function that is similar to setParameters.m which calls in your initial conditions. It is important that the order of these initial conditions match the order of your equations in step 2.

5. Adjust plotModelSimulation.m (if you want to simulate your model)
If you want to simulate your model, running this script will do so. "% ###Step 1" is setting custom colors for plotting purposes.
"% ###Step 2" is where the timespan for your model simulation is created, you can adjust t0 to change the initial time, and adjust tfinal to change the final time.
"% ###Step 3" sets options for the ODE solver, we used tolerances of 1e-12. You can change this if you want to adjust the accuracy of your solution (i.e. in case you want a different convergence speed), but this is not necessary.
"% ###Step 4" take the solutions to each differential equation and separates them by variable, to plot the solutions of each equation. You may want to adjust the names (left hand side of equal signs) to match your own variable names.
Everything after "% ###Step 5" makes a tiled plot layout of each simulation. For each plot, you will want to adjust the name of the Y variable if you adjusted the names in step 5e. Also adjustable are the colors, as mentioned in step a. You may also want to change the xlabel, ylabel, and title of each plot. If you have a different number of equations than the example, you may want to adjust "tiledlayout(2,3)", which is set to give you 2 rows and 3 columns of a tiled layout.

6. Adjust calculateQOI.m
Start with "% ###Step 1", adjusting for a different time span you want your model to be solved over. You can adjust t0 for a different initial time, and tfinal for a different final time.
"% ###Step 2" is where the sensitivity analysis will put the model outputs from. Y(end, 1) for example, pulls the value of the 1st equation at the last time point of the simulation. This last time point will likely match your tfinal time, as long as the solver is able to solve the system of equations within that timespan. Therefore, you only need to adjust “1” to match which equation in your odefun.m represents the variable that is your quantity of interest. For example, if your second equation is the output you are interested in, you should adjust this to Y(end, 2).

7. Adjust sobolMain.m
SobolMain is the file you want to run if you want to conduct a sensitivity analysis on your model. 
"% ###Step 1" sets an output folder that will save all files upon running SobolMain. This is helpful to refer to later if you need to keep track of your sensitvitiy analysis runs. The sensitivity indices will also be saved in this folder. Each file associated with a run will contain a timestamp in the form of YYYYMMDD_HHMMSS.
"% ###Step 2" will adjust how the range for each parameter is created. Each parameter will have a range [parameter*lower_percentage, parameter*upper_percentage]; you can adjust these as you see fit.
This step also determines how many base samples you want the sensitivity analysis to have. Note that if you are calculating S1 and ST values, the total number of model evaluations will be base_samples*(number of parameters + 2) and this of course will require more time and memory as you increase base_samples.
Finally, this step determines which type of probability distribution you want each parameter to have. You can adjust this.
That is all you need to adjust. Running this file now will output a few figures. The first figure shows scatterplots of the relationships between each parameter’s values and the QOI values. This is currently commented out for speed purposes, but is called 
% 4. Scatter Plot of QOI vs selected parameters". The second figure is a grouped bar graph, in descending order, ordered by ST that shows the first order and total order sensitivity indices side by side, with the parameter names on the x axis. You are welcome to adjust this and you can find the plotting lines in "% ###Step 3".
The rest of the files don’t need adjusting in order to run the sensitivity analysis. Optionally in "% ###Step 4" is a line that you can uncomment to plot just a few of the influential parameters; we used this to plot the ones with essentially non-zero indices, which we manually counted as the top 16. This is the end of what is needed to run the sensitivity analysis, but the following instructions are if you want to overlay histograms of one parameter set vs another. We used this to plot the distribution of the QOI values outputted from letting all parameters vary vs the distribution of the QOI values outputted from only letting the top influential parameters vary. Also listed are instructions for using makeSensitivityTables.m and plotsVaryingInfluentialParams.m.

9. Adjust compareQOIDistributions.m
"% ###Step 1" should match the previous settings from SobolMain.m.
"% ###Step 2" should define the parameter sets you want to use to compare QOI distributions outputted from letting only some parameters vary. We used all parameters as the first set, the top six influential parameters (known beforehand from running sobolMain.m), and then the least 29 influential parameters. You can have more than three sets, or fewer. The title here will be used in the legend of the comparison plot later on. The parameter names of the ones you are interested in letting vary should match setParameters.m. "% ###Step 3" is optional and is just to show you where the plots are configured. You can change the colors, displaystyle, etc.

10. makeSensitivityTables.m
This file conveniently makes a table for each type of sensitivity index. For each table, the files (saved from running SobolMain.m) are loaded in and then renamed. Then each table is populated and then sorted according to descending order of the run with the highest number of base samples (manually set, automatic detection of highest sample used is not present at the moment). We used this for the supplementary material.

11. plotsVaryingInfluentialParams.m
This file plots one at a time effects of influential parameters and then collective effects from tweaking all influential parameters to be in favor of the QOI and then all influential parameters to be against the QOI, as seen in the supplementary figures of the paper. % "% ###Step 1" has colors that can be adjusted. "% ###Step 2" has simulation settings that can be adjusted, including initial time, t0, and final time, tfinal, as well as tolerances in options. "% ###Step 3" currently solves the model six times, each time by tweaking only the parameter listed as the first argument of solveWithModifiedParam. "% ###Step 4" plots these simulations showing individual effects of changing one parameter on the QOI. "% ###Step 5" takes all parameters (manually) from step 4 and tweaks them all at once using plotAllForAllAgainst.
