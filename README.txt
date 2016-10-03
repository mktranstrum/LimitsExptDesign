The scripts in this folder should be used in the following order:
1.  Compile daskr, and models
./Compile.sh
2.  Regenerate data from brown
RegenerateBrownData.py
data saved in folder browndata
3.  Fit Brown data with mass action model
FitBrownDatatoMA.py seed
where seed is the value that will seed the random number generator
Note that this script requires geodesiclm package from https://sourceforge.net/projects/geodesiclm/
Edit this file to control the strength of the regularizing "Prior(s)"
parameter values are saved in folder fits
A good fit will have a Cost (1/2 Chi-squared) of a little more than 800.
4.  Generate artificial data to Apgar expts using MA model with fit parameters
GenerateApgarData.py
Edit this file to determie which parameter values (from step 3) are loaded
data saved in folder apgardata
5.  Fit Apgar data with Michael Menten model
FitApgarDatatoMM.py seed
where seed is the value that will seed the random number generator
Edit this file to determine the strength of the regulariging "Prior" and where it is centered
Typical fits have a Cost (1/2 Chi-squared) of ~50,000.
