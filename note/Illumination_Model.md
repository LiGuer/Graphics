
* Illumination Model
  - Define  
    Rendering equation:
    $$L_o(\boldsymbol p, \boldsymbol \omega_o, \lambda, t) = L_e(\boldsymbol p, \boldsymbol \omega_o, \lambda, t) + \int_\Omega L_i(\boldsymbol p, \boldsymbol \omega_i, \lambda, t) \cdot f_r(\boldsymbol p, \boldsymbol \omega_i, \boldsymbol \omega_o, \lambda, t) \cdot (\boldsymbol \omega_i \cdot \boldsymbol n) \mathrm d \boldsymbol \omega_i$$

    |Symbol|Mean|
    |---|---|
    | $\boldsymbol p$ | location in space |
    | $\lambda$ | particular wavelength of light |
    | $t$ | time |
    | $\boldsymbol n$ |surface normal at $p$|
    | $\boldsymbol \omega_o$ |direction of the outgoing light|
    | $\boldsymbol \omega_i$ |negative direction of the Incident light|
    | $\Omega$ |unit hemisphere centered around $\boldsymbol n$ containing all possible values for $\boldsymbol \omega_i$ where $\boldsymbol \omega_i \cdot \boldsymbol n > 0$|
    | $L_o$ |Total Outgoing light|
    | $L_e$ |Emitted light|
    | $L_i$ |Incident light|
    | $f_r$ |Bidirectional Reflectance Distribution function that defines how light is reflected at surface $\boldsymbol p$ from direction $\boldsymbol \omega_i$ to $\boldsymbol \omega_o$|
    |||


  - include 
    * Phong Model
      $$I_p = k_a I_a + \sum_{l \in lights} k_d (\hat {\boldsymbol L}_{i,l} \cdot \hat {\boldsymbol n}) I_{d,l} + k_{s,l} \left(\hat {\boldsymbol L}_{r,l} \cdot \hat {\boldsymbol v}_{eye}\right)^\alpha I_{s,l}$$

      $$\hat {\boldsymbol L}_{r,l} = 2 (\hat {\boldsymbol L}_{i,l} \cdot \hat {\boldsymbol n}) \hat {\boldsymbol n} - \hat {\boldsymbol L}_{i,l}  \tag{Reflection law}$$

      |Symbol|Mean|
      |---|---|
      | $k_a$ | ambient reflection constant |
      | $k_d$ | diffuse reflection constant |
      | $k_s$ | specular reflection constant |
      | $\alpha$ | shininess constant for this material |
      | $\hat {\boldsymbol n}$ | surface normal |
      | $\hat {\boldsymbol L}_{i,l}$ | negative direction of Incident light $l$ |
      | $\hat {\boldsymbol L}_{r,l}$ | direction of Reflect light of $l$ |
      | $\hat {\boldsymbol v}_{eye}$ | direction towards viewer |
      |||
      