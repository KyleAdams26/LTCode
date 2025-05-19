This repository contains the code used in "Understanding Immune Dynamics in Liver Transplant Through Mathematical Modeling" by Bruner et. al. The following are instructions on how to use the MATLAB scripts to run your own sensitivity analysis and model simulations using your own model.

The following document is a walkthrough of how to use the MATLAB code to conduct a sensitivity analysis on your own model using the code from the Liver Transplant Rejection project.

1. Download all scripts
SobolMain.m
compareQOIDists.m
getInitialConditions.m
getSamplesDesiredDist.m
getSamplesSobol.m
modelEQ.m
odefun.m
parameters.m
qoi.m
subset.m
sobolAnalysis.m

2. Adjust odefun.m to your model.
Start by adjusting lines 15-20. The vector “c” represents your state variables, and you can adjust the length of this vector here if needed. If you model as x_1, … , x_n state variables, you will have n lines here of x_i = c(i) where you can name x_i something meaningful.
Next, see lines 28-45. Here we have added the “pathways” of our model. Our pathways here are subsections of our equations. These pathways are made up of the state variable names you made in 2a and the parameters from the required structure “p” from the inputs of odefun.m. This “p” structure will have the names of our parameters along with their nominal values, and will be created in a separate file. Each parameter here should follow the example shown, in the form of p.parameterName.
Next, see lines 49 to 54. These should represent the dynamics/equations of your model. It is important these are in the same order as the state variables declared in step 2a. If in step 2a you declared x1 = c(1), then in this step, dy(1) must represent the derivative of x1 with respect to your independent variable. The number of these equations should also match the number of your state variables in step 2a. These dynamics should only refer to the parameters you have named in step 2c, and the variables you have named in step 2b. 
Lastly, adjust line 56 to contain as many equations as you made in step 2d.

3. Adjust parameters.m
Here you should create a line for each parameter, where the name of the parameter matches the names you used in step 2c. For each line, set p.parameterName = parameterValue, as shown in the example.
This script will be used anytime you need to load parameters into another script.

4. Adjust getInitialConditions.m
Here is a function that is similar to parameters.m which calls in your initial conditions. It is important that the order of these initial conditions match the order of your equations in step 2.
Adjust modelEQ.m (if you want to simulate your model)
If you want to simulate your model, running this script will do so. Lines 3-8 have custom colors made for plotting purposes. You are welcome to adjust these.
Line 11 loads in your parameters, you do not need to adjust this.
Line 14 is where the timespan for your model simulation is created, you can adjust t0 to change the initial time, and adjust tfinal to change the final time.
Line 20 sets options for the ODE solver, we used tolerances of 1e-12. You can change this if you want to adjust the accuracy of your solution (i.e. in case you want a different convergence speed), but this is not necessary.
Lines 25-30 take the solutions to each differential equation and separates them by variable, to plot the solutions of each equation. You may want to adjust the names (left hand side of equal signs) to match your own variable names.
Everything after line 35 makes a tiled plot layout of each simulation. For each plot, you will want to adjust the name of the Y variable if you adjusted the names in step 5e. Also adjustable are the colors, as mentioned in step a. You may also want to change the xlabel, ylabel, and title of each plot. If you have a different number of equations than the example, you may want to adjust line 35, which is set to give you 2 rows and 3 columns of a tiled layout.

5. Adjust qoi.m
Line 7 may be adjusted for a different time span you want your model to be solved over. You can adjust t0 for a different initial time, and tfinal for a different final time.
Line 16 is where the sensitivity analysis will put the model outputs from. Y(end, 6) for example, pulls the value of the 6th equation at the last time point of the simulation. This last time point will likely match your tfinal time, as long as the solver is able to solve the system of equations within that timespan. Therefore, you only need to adjust “6” to match which equation in your odefun.m represents the variable that is your quantity of interest. For example, if your first equation is the output you are interested in, you should adjust this to Y(end, 1).

7. Adjust SobolMain.m
SobolMain is the file you want to run if you want to conduct a sensitivity analysis on your model. 
Lines 9-10 will adjust how the range for each parameter is created. Each parameter will have a range [parameter*lower_percentage, parameter*upper_percentage], you can adjust these as you see fit.
Line 11 determines how many base samples you want the sensitivity analysis to have. Note that if you are calculating S1 and ST values, the total number of model evaluations will be base_samples*(number of parameters + 2) and this of course will require more time and memory as you increase base_samples.
Line 12 determines which type of probability distribution you want each parameter to have. You can adjust this.
Running this file now will output a few figures. The first figure shows scatterplots of the relationships between each parameter’s values and the QOI values. The second figure is a grouped bar graph, in descending order, ordered by ST that shows the first order and total order sensitivity indices side by side, with the parameter names on the x axis. You are welcome to adjust this and you can find the plotting lines within 82-94.
The rest of the files don’t need adjusting in order to run the sensitivity analysis. But the following instructions are if you want to overlay histograms of one parameter set vs another. We used this to plot the distribution of the QOI values outputted from letting all parameters vary vs the distribution of the QOI values outputted from only letting the top influential parameters vary.

9. Adjust compareQOIDists.m
Lines 9-12 should match step 7.
Line 13 should have the names of the parameters you want to let vary. These names should match exactly the names you used for parameters in all other scripts (parameters.m, ode.m, etc)
Line 14 should have how many parameters you’re varying, for the legend of the histogram that will be outputted if you run this file.
Line 68 sets up a 2-sample KS test to see how similar the distributions are.
