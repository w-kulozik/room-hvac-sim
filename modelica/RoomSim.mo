within;
package RoomSim
  "Room model of roomsim written independently in Modelica, for cross-validation (D58)"

  package Units
    type Temperature = Real(unit = "degC");
    type HeatCapacity = Real(unit = "J/K");
    type ThermalConductance = Real(unit = "W/K");
    type HumidityRatio = Real(unit = "kg/kg");
  end Units;

  function saturationPressure "Magnus form over liquid water (Alduchov and Eskridge 1996)"
    input Real T "Temperature, degC";
    output Real p "Saturation vapour pressure, Pa";
  algorithm
    p := 610.94*exp(17.625*T/(T + 243.04));
  end saturationPressure;

  function vapourPressure "Water-vapour partial pressure of moist air"
    input Real W "Humidity ratio, kg/kg";
    output Real p_w "Pa";
  protected
    constant Real P = 101325 "Total pressure, Pa";
    constant Real eps = 0.622 "Ratio of molar masses";
  algorithm
    p_w := W*P/(eps + W);
  end vapourPressure;

  function humidityRatio "Humidity ratio from temperature and relative humidity"
    input Real T "degC";
    input Real RH "Fraction";
    output Real W "kg/kg";
  protected
    constant Real P = 101325;
    constant Real eps = 0.622;
    Real p_w;
  algorithm
    p_w := RH*saturationPressure(T);
    W := eps*p_w/(P - p_w);
  end humidityRatio;

  model SummerDay "Sinusoidal summer day, the default scenario of roomsim"
    parameter Real T_mean = 26 "degC";
    parameter Real T_amp = 6 "K";
    parameter Real phase_h = 9 "h";
    parameter Real RH_o = 0.55 "Outdoor RH at T_mean, fraction";
    parameter Real G_max = 500 "Peak irradiance on the window, W/m2";
    parameter Real sunrise_h = 6;
    parameter Real sunset_h = 18;
    parameter Real n_people = 2;
    parameter Real occ_start_h = 8;
    parameter Real occ_end_h = 18;
    parameter Real C_o = 420 "ppm";
    output Real T_o "Outdoor temperature, degC";
    output Real W_o "Outdoor humidity ratio, kg/kg";
    output Real G_sol "Irradiance on the window, W/m2";
    output Real n_occ "Occupants";
    output Real h "Clock hour";
  equation
    h = mod(time, 86400)/3600;
    T_o = T_mean + T_amp*sin(2*Modelica.Constants.pi*(time/86400 - phase_h/24));
    W_o = humidityRatio(T_mean, RH_o);
    G_sol = if h >= sunrise_h and h < sunset_h then
      G_max*sin(Modelica.Constants.pi*(h - sunrise_h)/(sunset_h - sunrise_h)) else 0;
    n_occ = if h >= occ_start_h and h < occ_end_h then n_people else 0;
  end SummerDay;

  model Room "Air node, thermal mass, moisture with surface buffering, CO2"
    import Modelica.Constants.pi;
    // primary parameters (the defaults of roomsim.params.RoomParams)
    parameter Real volume = 5*4*2.7 "m3";
    parameter Real rho_air = 1.2 "kg/m3";
    parameter Real cp_air = 1006 "J/(kg K)";
    parameter Real furniture_factor = 5;
    parameter Units.HeatCapacity C_m = 10e6;
    parameter Real h_i = 3 "W/(m2 K)";
    parameter Real A_int = 60 "m2";
    parameter Units.ThermalConductance H_ao = 4;
    parameter Units.ThermalConductance H_mo = 25;
    parameter Real ach = 0.5 "1/h";
    parameter Real A_win = 3 "m2";
    parameter Real g_win = 0.6;
    parameter Real f_sol_air = 0.3;
    parameter Real q_person_sens = 70 "W";
    parameter Real f_person_conv = 0.5;
    parameter Real g_vap_person = 0.05/3600 "kg/s";
    parameter Real co2_person = 0.0052e-3 "m3/s";
    parameter Real buffer_area = 60 "m2";
    parameter Real buffer_density = 800 "kg/m3";
    parameter Real buffer_depth = 0.003 "m";
    parameter Real buffer_slope = 0.01 "kg/kg";
    parameter Real buffer_permeability = 2.5e-11 "kg/(m s Pa)";
    parameter Real T_ref = 25 "degC";
    // derived parameters, written from their definitions
    final parameter Real m_air = rho_air*volume "kg";
    final parameter Units.HeatCapacity C_a = rho_air*cp_air*volume*furniture_factor;
    final parameter Units.ThermalConductance H_am = h_i*A_int;
    final parameter Real V_inf = volume*ach/3600 "m3/s";
    final parameter Units.ThermalConductance H_inf = rho_air*cp_air*V_inf;
    final parameter Real C_b = buffer_density*buffer_depth*buffer_area*buffer_slope
      "Water stored per unit of buffer RH, kg";
    final parameter Real beta = h_i/(rho_air*cp_air)/(461.5*(T_ref + 273.15))
      "Surface vapour-transfer coefficient (Lewis analogy), kg/(m2 s Pa)";
    final parameter Real G_b = buffer_area/(1/beta + 0.5*buffer_depth/buffer_permeability)
      "Air-buffer vapour conductance, kg/(s Pa)";
    // inputs
    input Real T_o, W_o, G_sol, n_occ, C_o;
    input Real Q_hvac "Sensible heat added by the HVAC, W";
    input Real m_lat "Water removed by the HVAC, kg/s";
    input Real V_vent "Ventilation flow, m3/s";
    // states
    Units.Temperature T_a "Air temperature";
    Units.Temperature T_m "Mass temperature";
    Units.HumidityRatio W "Room humidity ratio";
    Real phi_b "RH in equilibrium with the buffer layer";
    Real C "CO2, ppm";
    Real m_buf "Water from the buffer to the air, kg/s";
    Real Q_sol "W";
  equation
    Q_sol = g_win*A_win*G_sol;
    m_buf = G_b*(phi_b*saturationPressure(T_m) - vapourPressure(W));
    C_a*der(T_a) = H_am*(T_m - T_a) + (H_ao + H_inf + rho_air*cp_air*V_vent)*(T_o - T_a)
      + f_sol_air*Q_sol + f_person_conv*n_occ*q_person_sens + Q_hvac;
    C_m*der(T_m) = H_am*(T_a - T_m) + H_mo*(T_o - T_m) + (1 - f_sol_air)*Q_sol
      + (1 - f_person_conv)*n_occ*q_person_sens;
    m_air*der(W) = rho_air*(V_inf + V_vent)*(W_o - W) + n_occ*g_vap_person - m_lat + m_buf;
    C_b*der(phi_b) = -m_buf;
    volume*der(C) = (V_inf + V_vent)*(C_o - C) + 1e6*n_occ*co2_person;
  end Room;

  model FreeFloating "Case A: three summer days without HVAC or ventilation"
    parameter Real T0 = 26 "Initial air and mass temperature, degC";
    SummerDay day;
    Room room(T_o = day.T_o, W_o = day.W_o, G_sol = day.G_sol, n_occ = day.n_occ,
      C_o = day.C_o, Q_hvac = 0, m_lat = 0, V_vent = 0);
  initial equation
    room.T_a = T0;
    room.T_m = T0;
    room.W = day.W_o;
    room.phi_b = vapourPressure(room.W)/saturationPressure(room.T_m);
    room.C = 420;
  end FreeFloating;

  model Scheduled "Case B: ventilation and fixed cooling and dehumidification while occupied"
    parameter Real T0 = 26;
    parameter Real Q_on = -800 "Sensible heat added while occupied, W";
    parameter Real m_on = 1e-4 "Water removed while occupied, kg/s";
    parameter Real V_on = 0.0105 "Ventilation while occupied, m3/s";
    SummerDay day;
    Boolean occupied = day.h >= 8 and day.h < 18;
    Room room(T_o = day.T_o, W_o = day.W_o, G_sol = day.G_sol, n_occ = day.n_occ,
      C_o = day.C_o, Q_hvac = if occupied then Q_on else 0,
      m_lat = if occupied then m_on else 0, V_vent = if occupied then V_on else 0);
  initial equation
    room.T_a = T0;
    room.T_m = T0;
    room.W = day.W_o;
    room.phi_b = vapourPressure(room.W)/saturationPressure(room.T_m);
    room.C = 420;
  end Scheduled;
end RoomSim;
