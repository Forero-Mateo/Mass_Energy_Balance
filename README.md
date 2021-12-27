# Undergrad M&E Balance II: Balancing with Recycle Stream
Project consisted of a 3-person team. All MATLAB scripting was done by me.

## Project Overview
This was a project for the junior level course Mass and Energy Balance II to solve a system of linear equations iteratively to solve an underdetermined system (recycle stream). Originally the project scope asked to solve for all unknown values for 3 different temperatures using Excels "Goal Seek" but having taught myself MATLAB the summer leading to this class, I chose to complete this project on MATLAB and solved for all 150 possible answers. 

The main purpose of demonstrating this project is to present the self-taught coding style and the numerical methods utilized (matrix algebra, solving underdetermined systems, advanced calculus etc.) to solve problems.

### Code Summary
The script is meant to be run all at once, so the first lines clear all stored variables, assigns an initial guess for the recycle stream, and sets the reactor temperature to 250⁰C. Following that, arrays are created to allocate memory for the plots that will be reported in the end. Different design equations are created to calculate the energy balance.

A while loop is used that runs iterates through different recycle stream guesses until an acceptable margin of error is reached. Once the mass balance is found, the energy balance is solved by calculating all relevant Q integral equations and the values are stored in the plot arrays. The temperature then increments by 1 and the initial recycle stream is set to 0 again. This loop continues until all 150 temperature set points are found and 6 subplots are output to analyze how reactor temperature impacts conversion %, energy usage, and recycle stream mass flow.
