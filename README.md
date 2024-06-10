# Model Calibration Steps for Norouzi et al. (2022)

This short note describes the steps required to calibrate the newly developed model in Norouzi et al. (2022). This repository contains two main MATLAB codes as well as two `.mat` files containing the soil retention measurements and spectral reflectance data.

## Contents

- `Main_Lebeau_Konrad_Model.m`
- `Main_Radiative_Transfer_Model.m`
- `soil_retention_data.mat`
- `spectral_reflectance_data.mat`

## Calibration Steps

### 1. Lebeau and Konrad (2010) Model Calibration

1. **Run the Code:**
   - Execute the `Main_Lebeau_Konrad_Model.m` script. This code fits the Lebeau and Konrad (2010) model to the measured soil water retention curve (SWRC).
   
2. **Optimization:**
   - The code saves the optimized parameters in a `.mat` file.
   - Set the variable `flag_optimization` to `0` if you do not need to repeat the optimization process each time you run the code. It will automatically load the previously saved results.

### 2. Radiative Transfer Model Calibration

1. **Run the Code:**
   - Execute the `Main_Radiative_Transfer_Model.m` script. This code loads the SWRC parameters and water components from the previous code as well as the hyperspectral reflectance measurements.
   
2. **Calibration:**
   - The code calibrates the new radiative transfer model [Eq. (8)] in the shortwave infrared range (i.e., 1300-2500 nm) via optimization of optical properties (`c_a`, `c_c`, `p_a`, `p_c`).
   - Set the variable `flag_opt` to `0` if you do not need to repeat the optimization process each time you run the code. The code will automatically load the previously saved results.

3. **Linear Mixing Model:**
   - If a linear mixing model (i.e., `p_a = p_c = 1`) is required, set the variable `linear_model_flag` to `1`.
   - For any value other than `1`, the code fits the general nonlinear mixing model.

4. **Model Predictions:**
   - In lines 17 and 18 of the script, you can specify the wavelengths and degree of saturations (`θ/θ_s`) at which model predictions are plotted.

## References

- Norouzi, S., Sadeghi, M., Tuller, M., Liaghat, A., Jones, S. B., & Ebrahimian, H. (2022). A novel physical-empirical model linking shortwave infrared reflectance and soil water retention. Journal of Hydrology, 614, 128653.
- Lebeau, M., & Konrad, J. M. (2010). A new capillary and thin film flow model for predicting the hydraulic conductivity of unsaturated porous media. Water Resources Research, 46(12).

## Contact

For any questions or issues, please contact sarem.nrz@gmail.com
