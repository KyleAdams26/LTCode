This repository contains the code used in "Understanding Immune Dynamics in Liver Transplant Through Mathematical Modeling" by Bruner et. al. in the Bulletin of Mathematical Biology, 2025. The instructions below explain how to modify the MATLAB scripts to run model simulations and sensitivity analysis for your own model.

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

2. Edit setParameters.m
In this file, replace the parameter names and their nominal values with those of your own model, one parameter per line. 
This script will be called anytime you need to load parameters into another script.

3. Edit odefun.m to your model.
a) In odefun.m, start with "% ###Step 1". The vector "c" represents your state variables. If your model has x_1, ..., x_n state variables, you will have n lines here of x_i = c(i) where you can name x_i something meaningful. You can have any number of state variables.
b) In "% $$$Step 2", the "pathways" of the model are defined. These can be used to simplify the equations that follow. The pathways are functions of the state variables and the parameters. Each parameter used in a pathway equation should follow the example shown, in the form of p.parameterName.
c) The section labeled "% ###Step 3" is where you write the equations of your model. These must be in the same order as the state variables declared in Step 1. If in Step 1 you declared x1 = c(1), then in Step 3, dy(1) must represent the derivative of x1 with respect to your independent variable. The number of equations should also match the number of your state variables in Step 1, and should only use parameters defined in setParameters.m and variables from Step 1.
d) Lastly, edit dcdt to contain as many equations as you made in Step 3.

4. Edit setInitialConditions.m
Edit this function to show your initial conditions. It is important that the order of these initial conditions matches the order of the equations in Step 3 of odefun.m .

5. Edit plotModelSimulation.m
If you want to simulate your model, edit and run this script.
a) "% ###Step 1" sets custom colors for each of the variables.
b) "% ###Step 2" sets the timespan for your model simulation. Edit t0 to change the initial time, and edit tfinal to change the final time.
c) "% ###Step 3" sets options for the ODE solver. We used tolerances of 1e-12. You can change this if you want to adjust the accuracy of your solution (e.g., if you want a different convergence speed), but this is not necessary.
d) "% ###Step 4" stores the solutions to each differential equation to be able to plot them later. Edit these variable names and numbers to match your model. 
e) Everything after "% ###Step 5" makes a tiled layout of separate plots of the simulated variables. For each plot, edit the name of the variable (after "Tf") to match the names in Step 4. Edit the colors (after 'Color') to match the color names from Step 1. Edit the xlabel, ylabel, and title of each plot to be appropriate for your model. If you have a different number of equations than the example, you can adjust "tiledlayout(2,3)", which is set to give you 2 rows and 3 columns of a tiled layout.

7. Edit calculateQOI.m
a) In "% ###Step 1", adjust t0 and tfinal for the time span you want your model to be solved over. 
b) "% ###Step 2" is where the quantity of interest (QOI) is saved after running the model at a set p of parameter values. Y(end, 1) pulls the value of the 1st equation at the last time point of the simulation. This last time point will likely match your tfinal time, as long as the solver is able to solve the system of equations within that timespan. You should adjust “1” to match the number of the equation in your odefun.m representing the quantity of interest. For example, if your second equation is the output you are interested in, you should change this to Y(end, 2).

8. Edit sobolMain.m
SobolMain runs Sobol sensitivity analysis on your model. 
a) "% ###Step 1" creates an output folder that will save all files when running SobolMain. This is helpful to refer to later if you need to keep track of your sensitvitiy analysis runs. The sensitivity indices will also be saved in this folder. Each file associated with a run will contain a timestamp in the form of YYYYMMDD_HHMMSS.
b) "% ###Step 2" sets the upper and lower bounds used to sample each parameter. Each parameter will have a range [parameter*lower_percentage, parameter*upper_percentage]; you can adjust these as needed.
This step also determines how many base samples to use for the sensitivity analysis. If you are calculating S1 and ST values, the total number of model evaluations will be base_samples*(number of parameters + 2) and this of course will require more time and memory as you increase base_samples. The code has base_samples set to a very low number (100) that can be used for initial testing of the code. 
Finally, this step specifies the type of probability distribution you want each parameter to have. You can adjust this.
Running this file now will output a few figures. The first figure will show scatterplots of the relationships between each parameter’s values and the QOI values. This is currently commented out for speed purposes, but is called "% 4. Scatter Plot of QOI vs selected parameters". The second figure is a grouped bar graph, in descending order, ordered by ST that shows the first order and total order sensitivity indices side by side, with the parameter names on the x axis. You can adjust this and you can find the plotting lines in "% ###Step 3".
d) The rest of the lines don’t need adjusting in order to run the sensitivity analysis. Optionally in "% ###Step 4" is a line that you can uncomment to plot just a few of the influential parameters; we used this to plot the ones with essentially non-zero indices, which for us was the top 16.

The following files create optional additional plots and tables.

10. Edit compareQOIDistributions.m
This file plots histograms of three sets of QOI values: the first set of QOI values is obtained by sampling and letting all parameters vary; the second only lets the top influential parameters vary; and the third one only lets the bottom least-influential parameters vary. These last two must be manually defined. 
a) "% ###Step 1" should match the previous settings from SobolMain.m.
b) In "% ###Step 2", list the three parameter sets whose QOI histograms will be plotted. We used all parameters as the first set, the top six influential parameters (determined from running sobolMain.m), and then the 29 least-influential parameters. You can have more than three sets, or fewer. The title here will be used in the legend of the comparison plot later on. The parameter names you list should match the versions of the names used setParameters.m.
c) Editing "% ###Step 3" is optional, but allows for changes to colors, displaystyle, etc.

12. Edit makeSensitivityTables.m
This file makes a table of values for each type of sensitivity index. For each table, the files (saved from running SobolMain.m) are loaded and then renamed. These file names must be manually entered in Step 1 (for Sobol total indices) and Step 2 (for Sobol first-order indices). Then each table is populated and then sorted according to descending order of the run with the highest number of base samples (this needs to be manually set, automatic detection of highest sample used is not present at the moment). We used this for the supplementary material.

13. Edit plotsVaryingInfluentialParams.m
This file plots the effects of one-at-a-time changes in influential parameters, as well as the results from changing all influential parameters to be in favor of the QOI and then all influential parameters to be against the QOI, as seen in the supplementary figures of the paper.
a) "% ###Step 1" has colors that can be adjusted. Brown is used because it is the color of our variable of interest, the L variable.
b) "% ###Step 2" has settings that can be adjusted, including initial time (t0) and final time (tfinal), as well as tolerances in "options".
c) "% ###Step 3" currently solves the model six times, each time by changing only the parameter listed as the first argument of solveWithModifiedParam (must be specified).
d) "% ###Step 4" plots these simulations showing individual effects of changing one parameter on the QOI. Parameter names must be specified.
e) "% ###Step 5" takes all parameters (specified manually) from Step 4 and changes them all at once using plotAllForAllAgainst.
