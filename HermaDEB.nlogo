; HermaDEB file entending the DEB-IBM Model of Ben Martin (btmarti25@gmail.com) and Elke Zimmer published in Martin, B. T., Zimmer, E. I., Grimm, V., & Jager, T. (2012). Dynamic Energy Budget theory meets individual‐based modelling: a generic and accessible implementation. Methods in Ecology and Evolution, 3(2), 445-449.

; Author: Dorra Louati : implementation of the standard DEB equations in an IBM with an extension to the hermaphrodite species (groupers)


; ==========================================================================================================================================
; ========================== DEFINITION OF PARAMETERS AND STATE VARIABLES ==================================================================
; ==========================================================================================================================================

; global parameters: are accessible for patches and turtles
globals[
  U_E^0    ; t L^2, initial reserves of the embryos at the start of the simulation
  f        ; - , scaled functional response
  L_0      ; m, initial structural volume
  month    ;
  year     ;
  i;
  j;
  y;
  m;
  flag ;
  filename ;

  ;cont;
]
; ------------------------------------------------------------------------------------------------------------------------------------------
; parameters for the environment: here only prey density

patches-own[
  X        ; # / m^2, prey density
  d_X      ; change of prey density in time
]
; ------------------------------------------------------------------------------------------------------------------------------------------

; definition of parameters for the individuals:
; the notation follows the DEBtool-notation as far as possible
; deviation: rates are indicated with "_rate" rather than a dot
; each individual(turtle) in the model has the following parameters
turtles-own[
  ; - - - - - - - - - - - - - - - STATE VARIABLES - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  L                ; m, structural length
  dL               ; change of structural length in time
  U_H              ; t L^2, scaled maturity
  dU_H             ; change of scaled maturity in time
  U_E              ; t L^2, scaled reserves
  dU_E             ; change of scaled reserves in time
  e_scaled         ; - , scaled reserves per unit of structure
  U_R              ; t L^2, scaled energy in reproduction buffer (not standard DEB)
  dU_R             ; change of energy in reproduction buffer (reproduction rate)
  sex              ; sexual role add in HermaDEB
  sex-change-list  ; sex-change historical list each item is
  offspring-list   ; off spring  historical list each item is
  coupled?         ; at each timestep this variable can change HermaDEB
  partner          ;
  offspring-count  ; with this parameter, the reproduction rate per turtle is shown on the interface
  annual-offspring-count
  inherited-trait  ; to be done later
  father ;
  mother ;
  potentiel-partner;
  sexually-active?
  age
  age-puberty
  l-puberty
  puberty-flag
  month-of-birth;
  year-of-birth;
  flag-sexchage
  first-change-age
  first-change-l
  l-group
  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ; - - - - - - - - - - - - - - - FLUXES (used by several submodels) - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  S_A         ; assimilation flux
  S_C         ; mobilisation flux

  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ; - - - - - - - - - - - - - - - EMBRYO (we use different state variable to not affect the state variable of the mother) - - - - - - - - --
  e_scaled_embryo
  e_ref
  U_E_embryo
  S_C_embryo
  U_H_embryo
  L_embryo
  dU_E_embryo
  dU_H_embryo
  dL_embryo
  ; parameters used to calculate the costs for an egg / initial reserves
  lower-bound ; lower boundary for shooting method
  upper-bound ; upper boundary for shooting method
  estimation  ; estimated value for the costs for an egg / initial reserve
  lay-egg?    ; parameter needed to hand over if an egg can be laid
  sim         ; this keeps track of how many times the calc-egg-size loop is run

  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ; - - - - - - - - - - - - - - - STANDARD DEB PARAMETERS (with dimension and name) - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  g           ; - , energy investment ratio
  v_rate      ; m /t , energy conductance (velocity)
  kap         ; - , allocation fraction to soma
  kap_R       ; - , reproduction efficiency
  k_M_rate    ; 1/t, somatic maintenance rate coefficient
  k_J_rate    ; 1/t, maturity maintenance rate coefficient
  U_H^b       ; t L^2, scaled maturity at birth
  U_H^p       ; t L^2, scaled maturity at puberty
  ; parameter that is used to randomize the input parameters
  scatter-multiplier

  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ; - - - - - - - - - - - - - - - PREY DYNAMICS (only relevant if prey-dynamics are set to logistic) - - - - - - - - - - - - - - - - - - -

  J_XAm_rate  ; # / (m^2 t), surface-area-specific maximum ingestion rate
  K           ; # / m^2, (half) saturation coefficient

  ; - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ; - - - - - - - - - - - - - - - AGEING -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  q_acceleration  ; - , ageing acceleration
  dq_acceleration ; change of ageing acceleration in time
  h_rate          ; - , hazard rate
  dh_rate         ; change of hazard rate in time
  age-day         ; each turtle has a random whole number between 0 and timestep if the mod of ticks = the age day of a turtle is will check to see if it dies
                  ; based on the ageing submodel. This is because mortality probabilities are per day, and timesteps are smaller
]

; ==========================================================================================================================================
; ========================== SETUP PROCEDURE: SETTING INITIAL CONDITIONS ===================================================================
; ==========================================================================================================================================

to setup
 ;; (for this model to work with NetLogo's new plotting features,
  ;; __clear-all-and-reset-ticks should be replaced with clear-all at
  ;; the beginning of your setup procedure and reset-ticks at the end
  ;; of the procedure.)
  __clear-all-and-reset-ticks
  ;file-close-all
 ;set filename  word replace-item 2 substring date-and-time 0 8 "h"  ".txt"
 ;set filename  word date-and-time  ".txt"
 set filename word "a"  ".txt"
  file-open  filename

 ;export-data-names-file
 export-short-data-names
 if add_my_pet? = "on"
 [convert-parameters]
 set i 0
 set j 0
 set y 0
 set m 0
 set month 1
 set year  0
 set L_0 .01  ; set initial length to some very small value (embryos start off as nearly all reserves)
 set flag 0


 crt 100; 100 turtles are created in the beginning
 ask  turtles  [
  individual-variability  ; first their individual variability in the parameter is set
  calc-embryo-reserve-investment     ; then the initial energy is calculated for each

 ]

ask patches [set X X_max] ; set initial value of energy per patch to X_max

end

; ==========================================================================================================================================
; ========================== GO PROCEDURE: RUNNING THE MODEL ===============================================================================
; ==========================================================================================================================================
; the go statement below is the order in which all procedures are run each timestep

to go


     if count turtles with [ U_H > U_H^p ] = 10 and flag = 0  [set flag 1
     set y year
     set m month
    ; show "j ai atteinds 20 adulte a" show month show year
    ]
   ask turtles
  [
    calc-dU_E                       ; first all individuals calculate the change in their state variables based on the current conditions
    calc-dU_H
    calc-dU_R
    calc-dL
    calc-age
  ]
  if aging = "on"                   ; if the ageing submodel is turned on, the change in damage inducing compound and damage are calculated
  [
    ask turtles
    [
      calc-dq_acceleration
      calc-dh_rate
      calc-age-L-puberty
    ]
  ]




if month = 3 [ask turtles with [U_H > U_H^p ] [
  randomlymove

  ]]


ask turtles with [U_H > U_H^p and (sex = "nosex" or sex = "female")] [

choose-sexual-role-by-firstchange-size
]


if month = one-of [6 7]
[

ask turtles  [ calc-lay-eggs  ]

ask turtles  with [sex = "female"]
 [
  if lay-egg? = 1 [
    find-partner-and-reproduce-female]
 ]


]

if month = 8 [              ask turtles with [sex = "female"] [set coupled? false ]]

ask turtles with [sex = "female"]
[move]



ask turtles with [U_H > U_H^p ]
[
 find-partner-neighbors
]
 calc_time
  update
  tick
  do-plots                          ; then the plots are updated

   if count turtles = 0  or year = 400 [  file-print mean [first-change-l] of turtles with [U_H > U_H^P]
     ifelse count turtles with [sex = "male"] > 0
     [file-print  Count turtles with [sex = "female"] / count turtles with [sex = "male"]]
     [file-print "zero male"]

;export-output filename
file-close
stop

     ]
end


;===========================================================================================================================================
;===========================================================================================================================================
;==============================Tools========================================================================================================



to calc_time
set i (ticks / timestep)
if (i mod 4) = 0
[
set month month + 1
;set i 0
]

if month = 13
[
;output-print year
ask turtles [
 export-data-output
 ;export-data-file
export-data-length
set annual-offspring-count 0
  ]



if year > 7 and count turtles with [ U_H > U_H^p] > 50 [ask n-of abs( count turtles with [ U_H > U_H^p] - 50)  turtles with [ U_H > U_H^p] [die]]
if year > 7 and count turtles with [ U_H > U_H^b and U_H < U_H^p] > 50 [ask n-of abs( count turtles with [ U_H > U_H^b and U_H < U_H^p]- 50)  turtles with [ U_H > U_H^b and U_H < U_H^p] [die]]


set year year + 1
set month 1
]
end


; ==========================================================================================================================================
; ========================== SUBMODELS =====================================================================================================
; ==========================================================================================================================================

; ---------------- conversion of parameters: from add_my_pet to standard DEB parameters ----------------------------------------------------

to convert-parameters
  let p_am p_m * zoom / kap_int
  set U_H^b_int E_H^b / p_am
  set U_H^p_int E_H^p / p_am
  set k_M_rate_int p_m / E_G
  set g_int (E_G * v_rate_int / p_am) / kap_int
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ------------------------ INDIVIDUAL VARIABILITY ------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------

to individual-variability
  ; individuals vary in their DEB paramters on a normal distribution with a mean on the input paramater and a coefficent of variation equal to the cv
  ; set cv to 0 for no variation
  set scatter-multiplier e ^ (random-normal 0 cv)
  set J_XAm_rate   J_XAm_rate_int * scatter-multiplier
  set g g_int / scatter-multiplier
  set U_H^b U_H^b_int / scatter-multiplier ;
  set U_H^p U_H^p_int / scatter-multiplier ;

  set v_rate v_rate_int
  set kap kap_int
  set kap_R kap_R_int
  set k_M_rate k_M_rate_int
  set k_J_rate k_J_rate_int
  set K  J_XAm_rate /   F_m
  set coupled? false
  set sex "nosex"
  set partner nobody
  set offspring-count 0
  ;set offspring-list
  ;set sex-change-list
  set annual-offspring-count 0
  set father nobody
  set mother nobody
  set potentiel-partner nobody
  set sexually-active? 0


  set shape "fish"
  set color 5
  set month-of-birth month
  set year-of-birth  year
  set age (list 0 0)
  set puberty-flag 0
  set age-puberty 0
  set l-puberty 0
  set flag-sexchage 0
  set first-change-age 0
  set first-change-l random-float 1.5 + Initial-length
  set l-group 0



end


; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- RESERVE DYNAMICS -------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; change in reserves: determined by the difference between assimilation (S_A) and mobilization (S_C) fluxes
; when food-dynamics are constant f = the value of f_scaled set in the user interface
; if food is set to  "logistic" f depends on prey density and the half-saturation coefficient (K)
; for embryos f = 0 because they do not feed exogenously

to calc-dU_E

  if food-dynamics = "constant"
  [ ifelse U_H <= U_H^b
    [set f 0]
   ; [set f f_scaled]
   [ifelse flag = 1 [set f X_moy] [set f f_scaled]]
  ]
  if food-dynamics = "logistic"
  [ ifelse U_H <= U_H^b
    [set f 0]
    [set f X / (K + X)]
  ]


  if food-dynamics = "fluctuation" or  food-dynamics = "step"
  [ ifelse U_H <= U_H^b
    [set f 0]
    [ifelse flag = 1 [set f X / (K + X)]  [set f f_scaled] ]
  ]
  set e_scaled v_rate * (U_E / L ^ 3)
  set S_C L ^ 2 * (g * e_scaled / (g + e_scaled)) * (1 + (L / (g * (V_rate / ( g * K_M_rate)))))

  set S_A  f * L ^ 2 ;

  set dU_E (S_A - S_C )
end
; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- MATURITY AND REPRODUCTION  ---------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; change in maturity is calculated (for immature individuals only)

to calc-dU_H

  ifelse U_H < U_H^p ; they only invest into maturity until they reach puberty
    [set dU_H ((1 - kap) * S_C - k_J_rate * U_H) ]
    [set dU_H 0]
end

; the following procedure calculates change in reprobuffer if mature
to calc-dU_R
  if U_H >= U_H^p
    [set dU_R  ((1 - kap) * S_C - k_J_rate * U_H^p) ]
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- DYNAMICS OF STRUCTURAL LENGHT-------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; the following procedure calculates change in structural length, if growth in negative the individual does not have enough energy to pay somatic maintenance and the starvation submodel is run
; where growth is set to 0 and individuals divirt enough energy from development (for juveniles) or reprodution (for adults) to pay maintenance costs
to calc-dL

  set dL   ((1 / 3) * (((V_rate /( g * L ^ 2 )) * S_C) - k_M_rate * L))

  if e_scaled < L / (V_rate / ( g * K_M_rate))  ; if growth is negative use starvation strategy 3 from the DEB book
    [
      set dl 0
      ifelse U_H < U_H^p
       [set dU_H (1 - kap) * e_scaled * L ^ 2 - K_J_rate * U_H^p - kap * L ^ 2 * ( L / (V_rate / ( g * K_M_rate)) - e_scaled)]
       [ set dU_R  (1 - kap) * e_scaled * L ^ 2 - K_J_rate * U_H^p - kap * L ^ 2 * ( L / (V_rate / ( g * K_M_rate)) - e_scaled)]
      set dU_E  S_A - e_scaled * L ^ 2
   ifelse U_H < U_H^p

 [  if dU_H < 0 [ ;create-juvenile
      die]  ]

      [if U_R < 0 [;create-juvenile
          die]]
    ]

end

;------------------------------------------------------------------------------------------------------------------------------------------
;---------- CHECK IF POSSIBLE TO LAY EGGS -------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------------------------------------------
; in the following, individuals determine if they have enough energy in their repro buffer to reproduce by creating an embryo with initial reserves set to the energy
; currently in their repro buffer * kap_R (conversion efficiancy of  reprobuffer to embryo) if the individual has enough energy to produce an offspring which will reach
; maturity and have a reserve density greater than the mothers when it hatches "lay-egg?" is set to one which will trigger the reproduction procedures "calc-egg-size" and "lay-eggs"
to calc-lay-eggs
  set L_embryo  L_0
  set U_E_embryo U_R * kap_R
  set U_H_embryo  0

  loop [
    set e_scaled_embryo v_rate * (U_E_embryo / L_embryo  ^ 3)
    set S_C_embryo L_embryo  ^ 2 * (g * e_scaled_embryo / (g + e_scaled_embryo)) * (1 + (L_embryo  / (g * (V_rate / ( g * K_M_rate)))))

    set dU_E_embryo  ( -1 * S_C_embryo )
    set dU_H_embryo  ((1 - kap) * S_C_embryo - k_J_rate * U_H_embryo )
    set dL_embryo  ((1 / 3) * (((V_rate /( g * L_embryo  ^ 2 )) * S_C_embryo) - k_M_rate * L_embryo ))

    set  U_E_embryo  U_E_embryo +  dU_E_embryo  / timestep
    set  U_H_embryo  U_H_embryo  +  dU_H_embryo   / timestep
    set  L_embryo    L_embryo  +  dL_embryo   / timestep

    if U_H_embryo  > U_H^b * 1 [ set lay-egg? 1 stop]
    if e_scaled_embryo < e_scaled [stop]
    ]
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ------------------------ INITIAL ENERGY --------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; calculate the initial energy of the first individuals using a bisection method

to calc-embryo-reserve-investment


  set lower-bound 0
  ifelse ticks = 0
  [set upper-bound 100]
  [set upper-bound U_R * kap_R]
  set sim 0

  loop[
    set sim sim + 1

    set estimation .5 * (lower-bound + upper-bound)
    set L_embryo  L_0
    set U_E_embryo estimation
    set U_H_embryo  0
    set e_scaled_embryo v_rate * (U_E_embryo / L_embryo  ^ 3)

    ifelse ticks = 0[set e_ref 1][set e_ref e_scaled]  ; e_ref now determines which e_scaled_embryo to calculate: 1 for ticks = 0 (in the setup procedure), e_scaled otherwise

    while [U_H_embryo  < U_H^b and e_scaled_embryo > e_ref ]
    ;     while [U_H_embryo  < U_H^b and e_scaled_embryo > 1 ] ; egg-size:  while [U_H_embryo  < U_H^b and e_scaled_embryo > e_scaled  ] ; I KEPT THIS LINE FOR NOW TO HAVE IT EASIER TO COMPARE
      [ set e_scaled_embryo v_rate * (U_E_embryo / L_embryo  ^ 3)
        set S_C_embryo L_embryo  ^ 2 * (g * e_scaled_embryo / (g + e_scaled_embryo)) * (1 + (L_embryo  / (g * (v_rate / ( g * k_M_rate)))))

        set dU_E_embryo  ( -1 * S_C_embryo )
        set dU_H_embryo  ((1 - kap) * S_C_embryo - k_J_rate * U_H_embryo  )
        set dL_embryo   ((1 / 3) * (((V_rate /( g * L_embryo  ^ 2 )) * S_C_embryo) - k_M_rate * L_embryo ))

        set  U_E_embryo  U_E_embryo +  dU_E_embryo    / (timestep )
        set  U_H_embryo   U_H_embryo  +  dU_H_embryo   / (timestep )
        set  L_embryo   L_embryo  +  dL_embryo    / (timestep )
      ]

    if e_scaled_embryo <  .05 +  e_ref and e_scaled_embryo > -.05 + e_ref and U_H_embryo  >= U_H^b  [

      ifelse ticks = 0 ;
      [set U_E^0 estimation
        set L L_0
        set U_E U_E^0
        set U_H 0
        set U_R 0
        set dU_R  0

        set age-day random timestep
        stop
      ][stop]]

    ifelse U_H_embryo  > U_H^b
      [ set upper-bound estimation ]
      [ set lower-bound estimation ]
    if sim > 100 [user-message ("Embryo submodel did not converge. Timestep may need to be smaller.") stop]
    ;if the timestep is too big relative to the speed of growth of species this will no converge
  ]
end


; added by Dorra 20-06-13 To direct set the initial values rather than use the dichotomous method

to calc-embryo-reserve-investment-modif

ifelse ticks = 0 ;
      [
       set U_E^0 0.01
        set L 0.036
        set U_E U_E^0
        set U_H 0
        set U_R 0
        set dU_R  0

        set age-day random timestep
        stop
      ][stop]

end
; end add




;-------------------------------------------------------------------------------------------------------------------------------------------
;--------- CHOOSE PARTNER SC1 ------------------------------------------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------------------------------------------------------------



to calc-age-L-puberty

  if  U_H > U_H^p and puberty-flag = 0
    [
      set age-puberty year
      set l-puberty L
      set puberty-flag 1
    ]

end
to choose-male-role
 let previous-sex sex
 set sex "male"
 set color 85
 if sex != previous-sex [set sex-change-list (list sex-change-list sex (list month year) l)]
end


to choose-sexual-role-bysize

let rlativ_f count turtles with [U_H > U_H^p and L <  [L] of myself]
let rlativ_m count turtles with [U_H > U_H^p and L >  [L] of myself]

if rlativ_m > rlativ_f  [set sex "female" set color 135]
if rlativ_f > rlativ_m  [set sex "male" set color 85]

end

to choose-sexual-role-by-threshold
   ifelse U_R > 10.98 [set sex "male" set color 85]  [set sex "female" set color 135]
end


to choose-sexual-role-by-sizethreshold
   ifelse L > L^* [set sex "male" set color 85]  [set sex "female" set color 135]
end

to choose-sexual-role-by-firstchange-size
   ifelse L >= first-change-l
   [set sex "male" set color 85
    set first-change-age last age]

   [set sex "female" set color 135]

end


to choose-sexual-role-byreserve

let previous-sex sex

ifelse potentiel-partner != nobody
[let rlativ_f count potentiel-partner with [ L <  [L] of myself]
let rlativ_m count potentiel-partner with [ L >  [L] of myself]

if rlativ_m > rlativ_f  [set sex "female" set color 135]
if rlativ_f > rlativ_m  [set sex "male" set color 85]
]
[  ifelse U_R > 10.98 [set sex "male" set color 85]  [set sex "female" set color 135]]
  ; set mylist lput 42 mylist
;set sex-change-list list [sex [month year] L] sex-change-list
if sex != previous-sex [set sex-change-list (list sex-change-list sex (list month year) l)]

end

to choose-sexual-role-bysexualdistribution

let previous-sex sex
find-partner-neighbors
ifelse potentiel-partner != nobody
[ifelse any? potentiel-partner with [L > [L] of myself and sex = "male"] ;< 0.5 * count potentiel-partner
  [set sex "female" set color 135]
  [set sex "male" set color 85]
]
[choose-sexual-role-bysize]
 if sex != previous-sex [set sex-change-list (list sex-change-list sex (list month year) l)]


 if previous-sex != "nosex" and sex != previous-sex and flag-sexchage = 0 [set first-change-age year
   set first-change-l l
   set flag-sexchage 1]


end




to find-partner-neighbors  ;; turtle procedure
  set potentiel-partner  other turtles with [U_H > U_H^p] in-radius vision
end

to find-partner-and-reproduce-male
  ; set partner one-of  turtles with [not coupled? and sex = "female" ] ;[not coupled? and sex != [sex] of myself and sex != "nosex"]
  set partner one-of  potentiel-partner with [not coupled? and sex = "female" ] ;[not coupled? and sex != [sex] of myself and sex != "nosex"]
  if partner != nobody
  [
        set coupled? true
        ask partner [ set coupled? true ]
        ask partner [ set partner myself ]

          ask partner [calc-embryo-reserve-investment-modif
                      lay-eggs
                      set partner nobody
                      ]
           set coupled? false
           set partner nobody



  ]

end

to find-partner-and-reproduce-female
  ; set partner one-of  turtles with [not coupled? and sex = "male" and L >  [L] of myself ] ;[not coupled? and sex != [sex] of myself and sex != "nosex"]
  if coupled? = false[
  find-partner-neighbors
  set partner one-of  potentiel-partner with [not coupled? and sex = "male" and L >  [L] of myself ] ;[not coupled? and sex != [sex] of myself and sex != "nosex"]
  if partner != nobody
  [
        set coupled? true
        set sexually-active? sexually-active? + 1
        ask partner [
          set coupled? true
          set partner myself
          set sexually-active? sexually-active? + 1
           ]

        calc-embryo-reserve-investment-modif
        lay-eggs-inheritance  ; individual inherit the average sex change length of their parents
        ask partner [ set coupled? false
                      set partner nobody]
        set partner nobody
  ]

  ]


end



;-------------------------------------------------------------------------------------------------------------------------------------------
;--------- LAY EGGS ------------------------------------------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------------------------------------------------------------
;the following procedure is run for mature individuals which have enough energy to reproduce
; they create 1 offspring and give it the following state variables and DEB parameters
;the initial reserves is set to the value determined by the bisection method in "calc_egg_size"

to lay-eggs

  let how-many-egg-laying int (U_R / estimation) * 0.01
  hatch how-many-egg-laying
    [
      ;the following code give offspring varibility in their DEB paramters on a normal distribution with a mean on the input paramater and a coefficent of variation equal to the cv
      ; set cv to 0 for no variation

      set scatter-multiplier e ^ (random-normal 0 cv)
      set J_XAm_rate   J_XAm_rate_int * scatter-multiplier
      set g g_int / scatter-multiplier
      set U_H^b    U_H^b_int / scatter-multiplier
      set U_H^p    U_H^p_int / scatter-multiplier

      set v_rate v_rate_int

      set kap kap_int
      set kap_R kap_R_int
      set k_M_rate k_M_rate_int
      set k_J_rate k_J_rate_int
      set  K J_XAm_rate /   F_m

      set L L_0
      set U_E estimation
      set U_H 0
      set U_R 0
      set dU_R  0
      set h_rate 0
      set dh_rate 0
      set q_acceleration 0
      set dq_acceleration 0
      set lay-egg? 0
      set age-day random timestep

      set coupled? false
      set sex "nosex"
      set partner nobody
      set offspring-count 0
      set mother myself
      set father [partner] of myself
    ]
  set offspring-count offspring-count + how-many-egg-laying
  ask partner [set offspring-count offspring-count + how-many-egg-laying * 0.01]
  set lay-egg? 0
  set U_R U_R - estimation * how-many-egg-laying
  ask partner [lay-sperm]
end

to lay-eggs-inheritance

  let eggnumber int ((U_R / estimation) * 0.000001 )
  let how-many-egg-laying eggnumber
  let half-how-many-egg-laying how-many-egg-laying / 2
  set offspring-count offspring-count + half-how-many-egg-laying
  set annual-offspring-count annual-offspring-count + half-how-many-egg-laying
  ;let how-many-egg-laying int (U_R / estimation) * 0.01
  let firstchange-l-mother first-change-l
  let firstchange-l-father  [first-change-l] of partner
  let mean-firstchange-l (firstchange-l-father + firstchange-l-mother) / 2
  let length-list list firstchange-l-father firstchange-l-mother
  hatch how-many-egg-laying
    [

      set scatter-multiplier e ^ (random-normal 0 cv)
      set J_XAm_rate   J_XAm_rate_int * scatter-multiplier
      set g g_int / scatter-multiplier
      set U_H^b    U_H^b_int / scatter-multiplier
      set U_H^p    U_H^p_int / scatter-multiplier

      set v_rate v_rate_int

      set kap kap_int
      set kap_R kap_R_int
      set k_M_rate k_M_rate_int
      set k_J_rate k_J_rate_int
      set  K J_XAm_rate /   F_m

      set L L_0
      set U_E estimation
      set U_H 0
      set U_R 0
      set dU_R  0
      set h_rate 0
      set dh_rate 0
      set q_acceleration 0
      set dq_acceleration 0
      set lay-egg? 0
      set puberty-flag 0
        set month-of-birth month
  set year-of-birth  year
      set age-day random timestep

      set coupled? false
      set sex "nosex"
      set partner nobody
      set offspring-count 0
      set mother myself
      set father [partner] of myself
      set first-change-age 0
     set first-change-l mean-firstchange-l
      ;
    ]
  set offspring-count offspring-count + half-how-many-egg-laying
  set lay-egg? 0
  set U_R U_R - estimation * how-many-egg-laying
   ask partner [lay-sperm
               set offspring-count offspring-count + half-how-many-egg-laying
               set annual-offspring-count annual-offspring-count + half-how-many-egg-laying

               ]


end


to lay-eggs-offspringonlycount
  let eggnumber int ((U_R / estimation) * 0.01 )
  let how-many-egg-laying eggnumber
  let half-how-many-egg-laying how-many-egg-laying / 2
  set offspring-count offspring-count + half-how-many-egg-laying
  set annual-offspring-count annual-offspring-count + half-how-many-egg-laying

  set lay-egg? 0
  set U_R U_R - estimation * how-many-egg-laying
  ask partner [lay-sperm
               set offspring-count offspring-count + half-how-many-egg-laying
               set annual-offspring-count annual-offspring-count + half-how-many-egg-laying

               ]
  end


to lay-sperm

  let how-many-sperm-laying int (U_R / estimation * sperm-clutch-cost)
   set U_R U_R - estimation * sperm-clutch-cost * how-many-sperm-laying

end



to create-juvenile

  hatch 1
    [
      ;the following code give offspring varibility in their DEB paramters on a normal distribution with a mean on the input paramater and a coefficent of variation equal to the cv
      ; set cv to 0 for no variation

      set scatter-multiplier e ^ (random-normal 0 cv)
      set J_XAm_rate   J_XAm_rate_int * scatter-multiplier
      set g g_int / scatter-multiplier
      set U_H^b    U_H^b_int / scatter-multiplier
      set U_H^p    U_H^p_int / scatter-multiplier

      set v_rate v_rate_int

      set kap kap_int
      set kap_R kap_R_int
      set k_M_rate k_M_rate_int
      set k_J_rate k_J_rate_int
      set  K J_XAm_rate /   F_m

      set L 0.48
      set U_E 6.39
      set U_H 4.44
      set U_R 0
      set dU_R  0
      set h_rate 4.8579145314179987E-4
      set dh_rate 1.5264792558146428E-6
      set q_acceleration 4.516462188831124E-6
      set dq_acceleration 2.9425183105396497E-9
      set lay-egg? 0
      set age-day random timestep
      set age  [7 8]


      set coupled? false
      set sex "nosex"
      set partner nobody
      set offspring-count 0
      set mother "god"
      set father "god" ;[partner] of myself
    ]


end

to show-juvenile-prop

 show L
 show U_E
 show U_H
 show U_R
 show dU_R
 show h_rate
 show dh_rate
 show q_acceleration
 show dq_acceleration
 show lay-egg?
 show age
; show age-day random timestep
end


; ------------------------------------------- ------------Update age -----------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
;

to calc-age

  let a year - year-of-birth
  let b month - month-of-birth
  set age (list b a)

end



; -------------------------------------------reintegrates the reserve-----------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
;For Males the unallocated reproductive energy reintegrates the reserve.

to U_R-reintegrates-e
  if  U_R > 0
  [set  U_E U_E + U_R]

end


; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- LOGISTIC PREY ----------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
 ;the following procedure calculates change in prey density this procedure is only run when prey dynamics are set to "logistic" in the user interface

to calc-d_X
   set d_X (r_X) * X * (1 - (X / K_X))   - sum [ S_A * J_XAm_rate   ] of turtles-here / volume
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- Food Fluctuation ----------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------


to calc-X-fluctuation
  let t ticks / timestep


  let A1 X_max - X_moy
  let A2 X_moy
  ifelse month = one-of [6 7 8 9]
  [
  set X A1 * sin( pi * t / 4 ) +  A2 + X_min
  ]

  [
  set X  A2 * -1 * sin( pi * (t - 4) / 8) +  A2 + X_min

  ]
end

to calc-X-step
  let t ticks / timestep


  let A1 X_max - X_moy
  let A2 X_moy
  ifelse month = one-of [6 7 8 9]
  [
  set X X_moy
  ]

  [
  set X X_min

  ]
end
; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- AGEING -----------------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
; the following procedure calculates the change in damage enducing compounds of an individual

to calc-dq_acceleration
  set dq_acceleration (q_acceleration * (L ^ 3 / (v_rate / ( g * k_M_rate)) ^ 3) * sG + H_a) * e_scaled * (( v_rate / L) - ((3 / L)*  dL)) - ((3 / L ) * dL) * q_acceleration
end

; the following procedure calculates the change in damage in the individual
to calc-dh_rate
  set dh_rate q_acceleration - ((3 / L) * dL) * h_rate
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- Export Data -----------------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------
 to export-data-output

   let status ""
   if U_H < U_H^b [set status "embryo"]
   if U_H > U_H^b and U_H < U_H^p [set status "juvenile"]
   if U_H > U_H^p [set status "adult"]
 ;  let l-group  ""
   if L < 0.25 [set l-group 0]
   if L > 0.25 and L < 0.35 [set l-group 1]
   if L > 0.35 and L < 0.45 [set l-group 2]
   if L > 0.45 and L < 0.55 [set l-group 3]
   if L > 0.55 and L < 0.65 [set l-group 4]
   if L > 0.65 and L < 0.75 [set l-group 5]
   if L > 0.75 and L < 0.85 [set l-group 6]
   if L > 0.85 and L < 0.95 [set l-group 7]
   if L > 0.95 and L < 1.05 [set l-group 8]
   if L > 1.05  [set l-group 9]
   output-write year
   output-write who
   output-write last age
   output-write age-puberty
   output-write l-puberty
  ; output-write ","
   output-write L
 ;  output-write ","
   output-write l-group
   output-write sex
 ;  output-write ","
   output-write status
 ;  output-write ","
   output-write annual-offspring-count
   output-write offspring-count
 ;  output-write ","
   output-write sexually-active?
   output-write first-change-age
   output-write first-change-l
   output-print ","

 end

  to export-data-file

   let status ""
   if U_H < U_H^b [set status "embryo"]
   if U_H > U_H^b and U_H < U_H^p [set status "juvenile"]
   if U_H > U_H^p [set status "adult"]
 ;  let l-group  ""
   if L < 0.25 [set l-group 0]
   if L > 0.25 and L < 0.35 [set l-group 1]
   if L > 0.35 and L < 0.45 [set l-group 2]
   if L > 0.45 and L < 0.55 [set l-group 3]
   if L > 0.55 and L < 0.65 [set l-group 4]
   if L > 0.65 and L < 0.75 [set l-group 5]
   if L > 0.75 and L < 0.85 [set l-group 6]
   if L > 0.85 and L < 0.95 [set l-group 7]
   if L > 0.95 and L < 1.05 [set l-group 8]
   if L > 1.05  [set l-group 9]
   file-write year
   file-write who
   file-write last age
   file-write age-puberty
   file-write l-puberty
  ; file-write ","
   file-write L
 ;  file-write ","
   file-write l-group
   file-write sex
 ;  file-write ","
   file-write status
 ;  file-write ","
   file-write annual-offspring-count
   file-write offspring-count
 ;  file-write ","
   file-write sexually-active?
   file-write first-change-age
   file-write first-change-l
   file-print ","

 end


 to export-data-length
   output-write year
   output-write who
   output-write last age
   output-write L
   output-write l-puberty
   output-write first-change-l
   output-print ","
 end

  to export-short-data-names
   output-write "year"
   output-write "who"
   output-write "age"
   output-write "L"
   output-write "l-puberty"
   output-write "first-change-l "
   output-print ","
  end
  to export-experiments_values
    output-write k_M_rate_int
   output-write g_int
   output-write  U_H^b_int
   output-write  U_H^p_int
   output-write h_a
  ; output-write ","
   output-write sG
 ;  output-write ","
   output-write v_rate_int
   output-write kap_int
 ;  output-write ","
   output-write kap_R_int
 ;  output-write ","
   output-write k_J_rate_int
   output-write vision
 ;  output-write ","
   output-write sperm-clutch-cost
    output-print ","


  end


 to export-data-names

    output-write "K_M_rate"
   output-write "g"
   output-write  "U_H^b"
   output-write  "U_H^p"
   output-write "ha"
  ; output-write ","
   output-write "SG"
 ;  output-write ","
   output-write "v"
   output-write "kap"
 ;  output-write ","
   output-write "kap_R"
 ;  output-write ","
   output-write "Kj"
   output-write "vision"
 ;  output-write ","
   output-write "sperm_cost"
    output-print ","
 export-experiments_values

   output-write "year"
   output-write "who"
   output-write  "age"
   output-write "age-puberty"
   output-write "l-puberty"
  ; output-write ","
   output-write "L"
 ;  output-write ","
   output-write "l-group"
   output-write "sex"
 ;  output-write ","
   output-write "status"
 ;  output-write ","
   output-write "annual-offspring-count"
   output-write "offspring-count"
 ;  output-write ","
   output-write "sexually-active?"
   output-write "first-change-age"
   output-write "first-change-l"
   output-print ","

 end

 to export-data-names-file

    file-write "K_M_rate"
   file-write "g"
   file-write  "U_H^b"
   file-write  "U_H^p"
   file-write "ha"
  ; file-write ","
   file-write "SG"
 ;  file-write ","
   file-write "v"
   file-write "kap"
 ;  file-write ","
   file-write "kap_R"
 ;  file-write ","
   file-write "Kj"
   file-write "vision"
 ;  file-write ","
   file-write "sperm_cost"
    file-print ","
 export-experiments-values-file

   file-write "year"
   file-write "who"
   file-write  "age"
   file-write "age-puberty"
   file-write "l-puberty"
  ; file-write ","
   file-write "L"
 ;  file-write ","
   file-write "l-group"
   file-write "sex"
 ;  file-write ","
   file-write "status"
 ;  file-write ","
   file-write "annual-offspring-count"
   file-write "offspring-count"
 ;  file-write ","
   file-write "sexually-active?"
   file-write "first-change-age"
   file-write "first-change-l"
   file-print ","

 end

  to export-experiments-values-file
    file-write k_M_rate_int
   file-write g_int
   file-write  U_H^b_int
   file-write  U_H^p_int
   file-write h_a
  ; file-write ","
   file-write sG
 ;  file-write ","
   file-write v_rate_int
   file-write kap_int
 ;  file-write ","
   file-write kap_R_int
 ;  file-write ","
   file-write k_J_rate_int
   file-write vision
 ;  file-write ","
   file-write sperm-clutch-cost
    file-print ","


  end

  to export-puberty-output
     if U_H > U_H^b and U_H < U_H^p
    [
    output-print "L and age for juvenil passage"
    output-write who
    output-write L
    output-write age
    ]
     if U_H > U_H^p
     [
    output-print "L and age for adult passage"
    output-write who
    output-write L
    output-write age
     ]
  end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- UPDATE -----------------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------

to update

  if year > 7
  [

    ask turtles with [U_H > U_H^b and L < 0.02]
  [
     if random-float 1 < background-mortality [ die]
  ]

  ]

  ask turtles
  [
    set U_E U_E + dU_E / timestep
    set U_H U_H + dU_H / timestep
    set U_R U_R + dU_R    / timestep
    set L L + dL    / timestep
    if U_H > U_H^b
    [ set q_acceleration q_acceleration + dq_acceleration  / timestep
      set h_rate h_rate + dh_rate  / timestep
    ]

   if aging = "on" [if ticks mod timestep = age-day [if random-float 1 < h_rate [die]] ] ;ageing related mortality
   if aging = "off" [if ticks mod timestep = age-day [if random-float 1 < background-mortality [

        die]] ]
 ]

  if food-dynamics = "fluctuation"[ ask patches [ calc-X-fluctuation]]
  if food-dynamics = "step"[ ask patches [ calc-X-step]]


end

to randomlymove  ;; turtle procedure randomly move
  setxy random-xcor random-ycor
end

to move

  rt random-float 360
  fd 1
end

; ------------------------------------------------------------------------------------------------------------------------------------------
; ----------------- PLOT -------------------------------------------------------------------------------------------------------------------
; ------------------------------------------------------------------------------------------------------------------------------------------

to do-plots
 set-current-plot "stage class density"
 set-current-plot-pen "embryo"
 set-plot-pen-interval 1 / timestep
    ifelse any? turtles with [U_H < U_H^b] [plot count turtles with [U_H < U_H^b]]
    [plot 0]
   set-current-plot-pen "juvenile"
   set-plot-pen-interval 1 / timestep
  ifelse any? turtles with [U_H > U_H^b and U_H < U_H^p] [plot count turtles with [U_H > U_H^b and U_H < U_H^p]]
  [plot 0]
  set-current-plot-pen "adult"
  set-plot-pen-interval 1 / timestep
  ifelse any? turtles with [U_H >= U_H^p] [plot count turtles with [U_H >= U_H^p]]
  [plot 0]

    set-current-plot-pen "male"
  set-plot-pen-interval 1 / timestep
  ifelse any? turtles with [U_H >= U_H^p] [plot count turtles with [U_H >= U_H^p and sex = "male"]]
  [plot 0]

    set-current-plot-pen "female"
  set-plot-pen-interval 1 / timestep
  ifelse any? turtles with [U_H >= U_H^p] [plot count turtles with [U_H >= U_H^p and sex = "female"]]
  [plot 0]


  set-current-plot "population density"
  set-plot-pen-interval 1 / timestep
 plot count turtles with [U_H > U_H^b]


end
@#$#@#$#@
GRAPHICS-WINDOW
1336
23
1581
261
11
11
9.0
1
14
1
1
1
0
1
1
1
-11
11
-11
11
1
1
1
ticks
30.0

SLIDER
279
575
391
608
f_scaled
f_scaled
0
1
1
.01
1
NIL
HORIZONTAL

BUTTON
42
70
108
103
go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
42
35
108
68
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
42
104
108
137
go-once
go
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
489
93
1311
294
stage class density
NIL
NIL
0.0
100.0
0.0
110.0
true
true
"" ""
PENS
"embryo" 1.0 0 -16777216 true "" ""
"juvenile" 1.0 0 -13345367 true "" ""
"adult" 1.0 0 -2674135 true "" ""
"male" 1.0 0 -817084 true "" ""
"female" 1.0 0 -5825686 true "" ""

SLIDER
110
35
243
68
timestep
timestep
0
1000
10
1
1
NIL
HORIZONTAL

MONITOR
264
10
321
55
Week
ticks / timestep
1
1
11

MONITOR
262
88
432
133
population density
count turtles
0
1
11

CHOOSER
279
530
458
575
food-dynamics
food-dynamics
"logistic" "constant" "fluctuation" "step"
2

INPUTBOX
39
467
120
527
v_rate_int
0.0178
1
0
Number

INPUTBOX
39
527
120
587
kap_int
0.8
1
0
Number

INPUTBOX
39
586
120
646
kap_R_int
0.95
1
0
Number

INPUTBOX
41
203
122
263
k_M_rate_int
0.0029
1
0
Number

INPUTBOX
40
646
120
706
k_J_rate_int
0.0025
1
0
Number

INPUTBOX
41
263
121
323
g_int
2.89
1
0
Number

INPUTBOX
42
323
121
383
U_H^b_int
5.22E-7
1
0
Number

INPUTBOX
42
383
121
443
U_H^p_int
4.4
1
0
Number

INPUTBOX
279
680
391
755
F_m
8.13
1
0
Number

INPUTBOX
391
575
458
635
r_X
0.5
1
0
Number

INPUTBOX
391
635
458
695
K_X
2
1
0
Number

INPUTBOX
391
695
458
755
volume
5
1
0
Number

INPUTBOX
278
608
391
680
J_XAm_rate_int
9.19E-4
1
0
Number

PLOT
489
294
1311
414
food density
NIL
NIL
0.0
2.0
0.0
2.0
true
false
"" ""
PENS
"pen-1" 4.0 0 -955883 true "" "plot  [X] of patch 0 0"

PLOT
490
414
1311
534
population density
NIL
NIL
0.0
100.0
0.0
2.0
true
false
"" ""
PENS
"> 2.6" 1.0 0 -2674135 true "" ""
"sex-ratio" 4.0 0 -11221820 true "" ""

CHOOSER
300
252
417
297
aging
aging
"on" "off"
0

INPUTBOX
300
297
417
357
h_a
1.0E-7
1
0
Number

INPUTBOX
300
357
417
417
sG
2.0E-4
1
0
Number

INPUTBOX
172
290
252
350
cv
0.1
1
0
Number

TEXTBOX
24
183
174
201
DEB-IBM parameters
11
0.0
1

TEXTBOX
304
509
454
527
feeding related parameters
11
0.0
1

TEXTBOX
318
229
468
247
ageing related parameters
11
0.0
1

INPUTBOX
160
465
255
525
p_m
0
1
0
Number

INPUTBOX
160
526
255
586
E_G
0
1
0
Number

INPUTBOX
159
706
253
766
zoom
2.61
1
0
Number

INPUTBOX
161
586
254
646
E_H^b
0
1
0
Number

INPUTBOX
159
646
254
706
E_H^p
0
1
0
Number

CHOOSER
160
420
255
465
add_my_pet?
add_my_pet?
"on" "off"
1

TEXTBOX
150
271
300
289
intraspecific variation
11
0.0
1

INPUTBOX
299
418
418
478
background-mortality
0
1
0
Number

MONITOR
578
38
703
83
NIL
mean [X] of patches
17
1
11

MONITOR
262
133
319
178
embryo
count turtles with [U_H < U_H^b]
17
1
11

MONITOR
318
133
375
178
juvenile
count turtles with [U_H > U_H^b and U_H < U_H^p]
17
1
11

MONITOR
374
133
431
178
Adult
count turtles with [ U_H > U_H^p]
17
1
11

MONITOR
321
10
378
55
Month
month
17
1
11

MONITOR
375
10
432
55
Year
year
17
1
11

MONITOR
318
178
375
223
Female
count turtles with [sex = \"female\"]
17
1
11

MONITOR
375
178
432
223
Male
count turtles with [sex = \"male\"]
17
1
11

SLIDER
110
73
239
106
sperm-clutch-cost
sperm-clutch-cost
0
1
0.1
0.1
1
NIL
HORIZONTAL

SLIDER
110
104
239
137
vision
vision
1
20
5
0.5
1
NIL
HORIZONTAL

BUTTON
41
136
108
169
Select and go
NIL
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
262
177
318
222
Sex-ratio
Count turtles with [sex = \"female\"] / count turtles with [sex = \"male\"]
17
1
11

INPUTBOX
483
583
638
643
X_max
1
1
0
Number

INPUTBOX
483
702
638
762
X_min
0.08
1
0
Number

INPUTBOX
483
642
640
702
X_moy
1
1
0
Number

TEXTBOX
510
563
660
581
Food Fluctuation
11
0.0
1

MONITOR
703
38
835
83
NIL
[X] of patch 0 0
17
1
11

MONITOR
487
38
579
83
NIL
count patches
17
1
11

OUTPUT
726
570
1221
756
11

INPUTBOX
172
350
252
410
L^*
0.7
1
0
Number

INPUTBOX
420
308
485
368
Initial-length
0.5
1
0
Number

TEXTBOX
279
62
416
80
Population statistics
12
0.0
1

TEXTBOX
605
14
755
32
Patches related\n
12
0.0
1

@#$#@#$#@
MODEL DESCRIPTION

A full model description, following the ODD protocol, is provided in "HermaDEB : An evolutionary IBM for energy allocation in hermaphrodites" which comes with this program.

## USER MANUAL

Species specific data is estimated using DEBTool and can be found here :

https://www.bio.vu.nl/thb/deb/deblab/add_my_pet/entries_web/Epinephelus_marginatus/Epinephelus_marginatus_res.html


## ZOOM

Depending on the resolution of the display you are using with your computer, you might not be able to see all elements of the Interface. In that case, please either use the scoll bars are the "Zoom" option in the menu.

## SPEED

You can speed up the program by deactivating the "view updates" option on the Interace tab.

## EXPORT DATA

To export model output:
- Use the file output primitives of NetLogo
- Right-click on the plots to export the data displayed to files
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
0
Rectangle -7500403 true true 151 225 180 285
Rectangle -7500403 true true 47 225 75 285
Rectangle -7500403 true true 15 75 210 225
Circle -7500403 true true 135 75 150
Circle -16777216 true false 165 76 116

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270

@#$#@#$#@
NetLogo 5.2.1
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>ticks = 1100</exitCondition>
    <metric>defects</metric>
    <enumeratedValueSet variable="timestep">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-sim">
      <value value="10"/>
      <value value="25"/>
      <value value="50"/>
      <value value="75"/>
      <value value="100"/>
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="bound-shift">
      <value value="0.05"/>
      <value value="0.1"/>
      <value value="0.2"/>
      <value value="0.3"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="zoom">
      <value value="0.2666"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_P_H">
      <value value="1.493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_k">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ref-t">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_U_b_H">
      <value value="9.993569821026685E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_B_H">
      <value value="0.0373"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_g">
      <value value="0.15003750937734434"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="-0.5"/>
      <value value="0"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_r">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_U_p_H">
      <value value="4.000107169649555E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_K_J_rate">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_kappa">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0.5"/>
      <value value="0.75"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="2800"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="T">
      <value value="20"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_kappa_r">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="0.00125"/>
      <value value="1.25E-5"/>
      <value value="1.25E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="arrhenius">
      <value value="6400"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Int_J_X_Am_rate">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="11200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_K_M_rate">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="int_v_rate">
      <value value="0.16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="100"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="sG">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="2.7698E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-Xmax" repetitions="1" runMetricsEveryStep="true">
    <setup>__clear-all-and-reset-ticks
 set i 0
 set j 0
 set y 0
 set m 0
 set month 1
 set year  0
 set L_0 .01  ; set initial length to some very small value (embryos start off as nearly all reserves) Dorra
 set flag 0

 crt 20; 10 turtles are created in the beginning
 ask  turtles  [
  individual-variability  ; first their individual variability in the parameter is set
  calc-embryo-reserve-investment     ; then the initial energy is calculated for each Dorra

 ]
ask patches [set X X_max]</setup>
    <go>go</go>
    <timeLimit steps="30000"/>
    <exitCondition>flag = 1 and  count turtles with [U_H &gt; U_H^P] &lt; 11</exitCondition>
    <metric>count turtles</metric>
    <metric>m</metric>
    <metric>y</metric>
    <metric>month</metric>
    <metric>year</metric>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="2.7698E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <steppedValueSet variable="X_max" first="0" step="0.1" last="2"/>
    <enumeratedValueSet variable="cv">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;Off&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="2.7698E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="1.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;Off&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentFinattestForcountLgroup" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>year</metric>
    <metric>month</metric>
    <metric>count turtles with [l-group = 0]</metric>
    <metric>count turtles with [l-group = 1]</metric>
    <metric>count turtles with [l-group = 2]</metric>
    <metric>count turtles with [l-group = 3]</metric>
    <metric>count turtles with [l-group = 4]</metric>
    <metric>count turtles with [l-group = 5]</metric>
    <metric>count turtles with [l-group = 6]</metric>
    <metric>count turtles with [l-group = 7]</metric>
    <metric>count turtles with [l-group = 8]</metric>
    <metric>count turtles with [l-group = 9]</metric>
    <metric>count turtles with [l-group = 0 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 1 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 2 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 3 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 4 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 5 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 6 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 7 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 8 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 9 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 0 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 1 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 2 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 3 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 4 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 5 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 6 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 7 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 8 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 9 and sex = "female"]</metric>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.4E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="5" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>year</metric>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <steppedValueSet variable="h_a" first="1.0E-6" step="1.0E-6" last="9.0E-6"/>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-annual_offspring-count-bysex" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>year</metric>
    <metric>month</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 1 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 2 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 3 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 4 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 5 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 6 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 7 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 8 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 9 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 1 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 2 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 3 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 4 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 5 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 6 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 7 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 8 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 9 and sex = "female"]</metric>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment-firstchagelength" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>year</metric>
    <metric>month</metric>
    <metric>[first-change-l] of turtles</metric>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentstepcount" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>year</metric>
    <metric>month</metric>
    <metric>count turtles with [l-group = 0]</metric>
    <metric>count turtles with [l-group = 1]</metric>
    <metric>count turtles with [l-group = 2]</metric>
    <metric>count turtles with [l-group = 3]</metric>
    <metric>count turtles with [l-group = 4]</metric>
    <metric>count turtles with [l-group = 5]</metric>
    <metric>count turtles with [l-group = 6]</metric>
    <metric>count turtles with [l-group = 7]</metric>
    <metric>count turtles with [l-group = 8]</metric>
    <metric>count turtles with [l-group = 9]</metric>
    <metric>count turtles with [l-group = 0 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 1 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 2 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 3 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 4 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 5 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 6 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 7 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 8 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 9 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 0 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 1 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 2 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 3 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 4 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 5 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 6 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 7 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 8 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 9 and sex = "female"]</metric>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;step&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentFluctuationcount" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>year</metric>
    <metric>month</metric>
    <metric>count turtles with [l-group = 0]</metric>
    <metric>count turtles with [l-group = 1]</metric>
    <metric>count turtles with [l-group = 2]</metric>
    <metric>count turtles with [l-group = 3]</metric>
    <metric>count turtles with [l-group = 4]</metric>
    <metric>count turtles with [l-group = 5]</metric>
    <metric>count turtles with [l-group = 6]</metric>
    <metric>count turtles with [l-group = 7]</metric>
    <metric>count turtles with [l-group = 8]</metric>
    <metric>count turtles with [l-group = 9]</metric>
    <metric>count turtles with [l-group = 0 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 1 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 2 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 3 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 4 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 5 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 6 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 7 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 8 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 9 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 0 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 1 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 2 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 3 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 4 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 5 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 6 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 7 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 8 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 9 and sex = "female"]</metric>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;fluctuation&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentoffspring_fluctuation_bysex" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>year</metric>
    <metric>month</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 1 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 2 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 3 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 4 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 5 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 6 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 7 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 8 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 9 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 1 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 2 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 3 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 4 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 5 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 6 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 7 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 8 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 9 and sex = "female"]</metric>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;fluctuation&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentoffspringcount-stepfunction" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>year</metric>
    <metric>month</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 1 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 2 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 3 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 4 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 5 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 6 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 7 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 8 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 9 and sex = "male"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 1 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 2 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 3 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 4 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 5 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 6 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 7 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 8 and sex = "female"]</metric>
    <metric>mean [annual-offspring-count] of turtles with [l-group = 9 and sex = "female"]</metric>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;step&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentforconstpopLocaldistsexchoicestep" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="100000"/>
    <metric>year</metric>
    <metric>month</metric>
    <metric>count turtles with [l-group = 0]</metric>
    <metric>count turtles with [l-group = 1]</metric>
    <metric>count turtles with [l-group = 2]</metric>
    <metric>count turtles with [l-group = 3]</metric>
    <metric>count turtles with [l-group = 4]</metric>
    <metric>count turtles with [l-group = 5]</metric>
    <metric>count turtles with [l-group = 6]</metric>
    <metric>count turtles with [l-group = 7]</metric>
    <metric>count turtles with [l-group = 8]</metric>
    <metric>count turtles with [l-group = 9]</metric>
    <metric>count turtles with [l-group = 0 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 1 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 2 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 3 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 4 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 5 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 6 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 7 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 8 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 9 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 0 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 1 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 2 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 3 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 4 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 5 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 6 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 7 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 8 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 9 and sex = "female"]</metric>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;step&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="L^*">
      <value value="0.7"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment2402constfoodstepfluctuation" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="100000"/>
    <metric>year</metric>
    <metric>month</metric>
    <metric>count turtles with [l-group = 0]</metric>
    <metric>count turtles with [l-group = 1]</metric>
    <metric>count turtles with [l-group = 2]</metric>
    <metric>count turtles with [l-group = 3]</metric>
    <metric>count turtles with [l-group = 4]</metric>
    <metric>count turtles with [l-group = 5]</metric>
    <metric>count turtles with [l-group = 6]</metric>
    <metric>count turtles with [l-group = 7]</metric>
    <metric>count turtles with [l-group = 8]</metric>
    <metric>count turtles with [l-group = 9]</metric>
    <metric>count turtles with [l-group = 0 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 1 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 2 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 3 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 4 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 5 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 6 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 7 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 8 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 9 and sex = "male"]</metric>
    <metric>count turtles with [l-group = 0 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 1 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 2 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 3 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 4 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 5 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 6 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 7 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 8 and sex = "female"]</metric>
    <metric>count turtles with [l-group = 9 and sex = "female"]</metric>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;step&quot;"/>
      <value value="&quot;fluctuation&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="L^*">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="0.12"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="L^*">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentLengthoffatheronly" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <exitCondition>year = 1200</exitCondition>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="L^*">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentMeanParents" repetitions="2" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <enumeratedValueSet variable="L^*">
      <value value="0.7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="J_XAm_rate_int">
      <value value="9.19E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="vision">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="food-dynamics">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^b">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^p_int">
      <value value="4.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_J_rate_int">
      <value value="0.0025"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_int">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="timestep">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="f_scaled">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sG">
      <value value="2.0E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="background-mortality">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="v_rate_int">
      <value value="0.0178"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="g_int">
      <value value="2.89"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cv">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="p_m">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_max">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="F_m">
      <value value="8.13"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_moy">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Reproduct">
      <value value="&quot;On&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="sperm-clutch-cost">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="h_a">
      <value value="1.0E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="X_min">
      <value value="0.08"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_G">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="add_my_pet?">
      <value value="&quot;off&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="E_H^p">
      <value value="0"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k_M_rate_int">
      <value value="0.0029"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="kap_R_int">
      <value value="0.95"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="U_H^b_int">
      <value value="5.22E-7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="volume">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="aging">
      <value value="&quot;on&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="r_X">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K_X">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="zoom">
      <value value="2.61"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180

@#$#@#$#@
0
@#$#@#$#@
