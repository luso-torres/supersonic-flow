# Numerical Simulation of a Supersonic Diffuser and Nozzle

This study aims to validate the optimal geometries derived from the analytical calculations of the diffuser and nozzle. 

Part I focuses on the diffuser, while Part II will assess the nozzles based on their inviscid and viscous behavior.

The results for the Mach Number and Velocity Vector (viscous case) are illustrated below. The full document contains information about the flow properties (Mach Number, Pressure, Temperature, Velocity and Density) along with graphs showing the change of these properties.

All simulations presented in this document were conducted using the tools provided by ANSYS, specifically the CFD application FLUENT.
For full details of the project, checkout the reports containing the methods.


## Geometry
The whole geometry is illustrated as follows:
![nozzle](https://github.com/user-attachments/assets/d011b4e2-c0a6-4e38-81c3-ada391aaeb54)

1. Oblique Shockwave on region 1.
2. Oblique Shockwave on region 2.
3. Normal Shockwave when entering the Diffuser 3.
4. Combustion Chamber (Heat addition plus increase in area) at 4.
5. Inlet of the of the Nozzle.
6. Transonic Region.
7. Exit of the Nozzle.

### Diffuser
The initial configuration of the Diffuser has three stages of shockwaves:
<img src="https://github.com/user-attachments/assets/d19b86e9-dd60-4e81-b9d0-d75ee2b4802a" alt="Diffuser" width="600" height="300">

The analytical results obtained in the first work are shown below:
| Parameter                | 0 (∞)  | 1     | 2     | 3     |
|--------------------------|---------|-------|-------|-------|
| Mach Number              | 3.5     | 2.55  | 1.65  | 0.66  |
| Static Pressure [kPa]    | 5.47    | 18.58 | 63.02 | 189.12|
| Temperature [K]          | 216.65  | 325.82| 489.97| 701.21|
| Velocity [km/s]          | 1.03    | 0.93  | 0.73  | 0.35  |
| Density [kg/m³]          | 0.09    | 0.20  | 0.45  | 0.94  |

Based on the isentropic properties of air in this configuration, the following values for the geometry then are derived and presented in the figure below:
| **Property**                           | **Value**  |
|---------------------------------------|------------|
| $\beta_1$ [deg]                       | 29.90      |
| $\beta_2$ [deg]                       | 43.10      |
| $\theta_1$ [deg]                      | 15.60      |
| $\theta_2$ [deg]                      | 20.57      |
| $\left(\frac{P_{0,e}}{P_{0,i}}\right)$ | 0.62       |
| $D$ [kN]                              | 1716.5     |

<img src="https://github.com/user-attachments/assets/0145468f-e8dd-4951-9245-f5481b8f3089" alt="Diffuser Geometry" width="800" height="400">

## Nozzles
In this work, two geometries of nozzles were investigated, being the first a conical-shaped nozzle, and in the second case a bell-shaped nozzle.


### Conical-shaped
<img src="https://github.com/user-attachments/assets/8889b24a-b99a-49c7-88c8-a751efc4476f" alt="Diffuser Mach Number Viscous" width="800" height="400">


### Bell-shaped

<img src="https://github.com/user-attachments/assets/116d25b7-0b75-4e40-8c78-cd9b459e44dc" alt="Diffuser Mach Number Viscous" width="800" height="400">



## Flow Behavior for the Viscous Case:
### Diffuser
#### Mach Number
<img src="https://github.com/user-attachments/assets/b448ccc7-e39b-4cc7-a283-6882b3b4e283" alt="Diffuser Mach Number Viscous" width="800" height="400">

#### Velocity Vector
<img src="https://github.com/user-attachments/assets/7edfaaec-d43f-4b52-b135-656345b34d1f" alt="Diffuser Velocity Viscous" width="800" height="400">

### Nozzles

#### Conical-Shape

#### Mach Number
<img src="https://github.com/user-attachments/assets/1f80da7b-3108-432c-a88c-c077b18c7535" alt="Conical Nozzle Mach Number Viscous" width="800" height="400">

#### Velocity Vector
<img src="https://github.com/user-attachments/assets/a4b91434-60ac-44fe-805c-26619fc3e8b8" alt="Conical Nozzle Velocity Vector Viscous" width="800" height="400">



#### Bell-Shaped
#### Mach Number
<img src="https://github.com/user-attachments/assets/81be9de5-916f-4b90-827a-0b6435f6c2b8" alt="Bell-Shaped Nozzle Mach Number Viscous" width="800" height="400">

#### Velocity Vector
<img src="https://github.com/user-attachments/assets/407b630d-ff53-419e-b9ff-c79fa63e95c4" alt="Bell-Shaped Nozzle Velocity Vector Viscous" width="800" height="400">






