;extensions[vid]

globals
[
  ;;; landscape design
  proportionCrops ; list of the different crop proportion
  orderCrops ; if they are distributed around rops it is on this order

  totalPatchesLandscape ; total number of patches in the landscape
  totalSNHPatches
  totalCropPatches ; total number of crop patches in the landscape
  totalCiPatches ; list of the total number of crop i in the landscape

  ;;; time
  date ; date during the year
  totalTicksSimulation ; total nb. of ticks = nb. years * length-season
  year ; year counter

  ;;; pest dynamic
  sensibilityFrames ; list of the different frames of sensibilities [[0 50] [60 120] ....]
  cropsInfectionrates; list of the differents infection rates
  nbOfCiToInfectPerTick; list of the differents numbers of patches that have to be infected each tick.
  theoreticalNbFinallyInfectedCrops
  infectionFrequency

  ;;; adult-predators dynamic
  detectionRadius

  densityDependenceMortalityParameter ; for pattern of mortality with density-dependence

  totalDeaths ; total nb. of adult-predators which die during a year
  totalMovements ; total nb. of foraging movements by adult-predators
  totalMovementsInCrops ; total nb. of foraging + overwintering movements by adult-predators in CROPS
  totalMovementsInSNH ; total nb. of foraging movements by adult-predators in SNH

  ;;; files for outputs
  fileName ; var. root file-name

  ; var. file names
  tickFileName
  endSeasonFileName
  spatialFileName
  movementsFileName

  ;;; indicators
  regulation

]
;;;;;;; files to include
__includes
[
  "methods.nls"
  "landscape-design.nls"
  "crops-pests-dynamics.nls"
  "adult-predators-dynamic.nls"
]
;;;;;;; patches and breed
patches-own[
  ;;; landscape
  landUse ; crop or SNH (1 or 0)
  landCover ; SNH and the different crops (0 1 2 3 4 5 ...)
  fixed ; used when crops have to be distributed around SNH

  cluster ; identity of the leader patch of my own cluster
  numCluster

  closestSNHDistance ; distance (nb of patches) to the closest SNH patch
  xClosestSNH ; pxcor of the closest SNH patch
  yClosestSNH ; pycor of the closest SNH patch

  ;;; status

  state ; state = sane, latent, attractive, with juveniles
  nextState

  predatorPresence ; is there an adult-predator or a juvenile on this patch?

  ;;; dynamic

  infectionDate ; date of pest arrival
  infectionRate
  sensibilityD1 ; beginning of their frame sensibility
  sensibilityD2 ; end of the frame sensibility

  lengthLatentPeriod ; how many ticks before being attractive?

  datePredatorArrival
  firstPredatorArrived ; who is the first predator on me?
  birthYearOfMyPredator ; to know if the predator that stopped on this crop was of year n or year n-1
  timeForCropLoss ; duration between pest and predator arrivals

  ;;; only for crop patches
  foraging-movement-distance-for-predator-before-arrival ; the distance that the predator has done, in both snh or crops, before arrival

  durationBeforeEggs ; how many ticks before the adult-predator on me lay eggs?
  nbTicksSinceAdultArrival ; nb. of tick since adult arrival

  lenghtJuvenileMaturation ; how many ticks before the juvenile on me converts itself in adult?
  nbTicksSinceJuvenileBirth ; nb. of ticks since juvenile birth

  ;;; outputs

  visitCounter
  totalInfectionCycles ; count nb. of times a patch is infected
  totalCurationCycles ; count cycles of SIS (for a single patch)
  yield
]

breed [adultPredators adultPredator]

adultPredators-own
[
  birthDate ; date the predator was born

  ;;;
  birthYear ; year the predator was born

  foragingValue ; variable that will indicate when the predator move (if foragingMoveFrequency != 1)

  ;;; date for starting foraging-movement
  dateFirstForagingMovement ; date that allows the adult-predator to forage

  ;;;
  distanceInPatches

  distanceForOverwintering ; distance travelled during overwintering movement

]

breed [juvenilePredators juvenilePredator]


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; main
to setup
  clear-all
  reset-ticks

  initializeParametersMethod
  landscapeDesignMethod

  ask patches [
   ifelse landCover > 0
   [
     set sensibilityD1 item 0 (item (landCover - 1) sensibilityFrames)
     set sensibilityD2 item 1 (item (landCover - 1) sensibilityFrames)
     set infectionRate item (landCover - 1) cropsInfectionRates
   ]
   [
     set sensibilityD1 0
     set sensibilityD2 0
   ]
  ]

  set totalTicksSimulation (totalYearsSimulation * lengthSeason)

  initializeAdultPredatorsMethod
  set date 0
  set year 1
  set-current-directory (word "/home/tess/Desktop/Data_Stage" folderPath)
  initializeFileNamesMethod

end

to setup2
  clear-all
  reset-ticks

  initializeParametersMethod

  if config = "haut gauche"
  [
    ask patches with [pxcor + pycor < 0][set landUse 1 set landCover 1 set fixed TRUE]
    ask patches with [pxcor + pycor >= 0][set landUse 1 set landCover 3 set fixed TRUE]
    ask patches with [pxcor - pycor >= 8][set landUse 1 set landCover 2 set fixed TRUE]
    ask patches with [pxcor < -5 and pycor > 6][set landUse 0 set landCover 0 set fixed TRUE]
  ]

  if config = "bande"
  [ ask patches with[pxcor >= -13] [set landUse 1 set landCover 1 set fixed TRUE]
    ask patches with [pxcor  >= -3 ][set landUse 1 set landCover 2 set fixed TRUE]
    ask patches with [pxcor >= 7][set landUse 1 set landCover 3 set fixed TRUE]
    ask patches with [pxcor < -13][set landUse 0 set landCover 0 set fixed TRUE]
  ]
  if config = "centre123"
  [
    let counter 1
    let rajout 0
    let objective precision (proportionSNHPatches / 100 * count patches) 0
    print((word "nb de patch souhaite " objective))
    let distance_centre 1

    ask patches [set landUse 1 set landCover (last orderCrops) set pcolor 55 - 10 * last orderCrops]
    ask patch 0 0 [set landUse 0 set landCover 0 set fixed TRUE]

    while [counter < objective] [
      set rajout rajout + 8
      ifelse objective - counter >= rajout [
        ask patches with [fixed != TRUE and abs pxcor <= distance_centre and abs pycor <= distance_centre][set landUse 0 set landCover 0 set pcolor 55 set fixed TRUE]
        set counter counter + rajout
        set distance_centre distance_centre + 1
      ][
        ask n-of (objective - counter) patches with [abs pxcor <= distance_centre and abs pycor <= distance_centre and fixed != TRUE][set landUse 0 set landCover 0 set pcolor 55 set fixed TRUE]
        set counter counter + (objective - counter)
      ]

    ]
    print((word "nombre de SNH:" count patches with [landCover = 0]))

    set objective precision ((1089 - objective) / 3) 0
    print((word "nb de patch par culture visé: " objective))
    (foreach butlast orderCrops
      [[x] ->
        print((word"-----crop " x " -----"))
        ask patches with [abs pxcor <= distance_centre and abs pycor <= distance_centre and fixed != TRUE][set landUse 1 set landCover x set pcolor 55 - 10 * x set fixed TRUE]
        set counter count patches with [landCover = x]
        print((word "pour completer " count patches with [landCover = x]))

        while [counter < objective][
          set rajout rajout + 8
          ifelse objective - counter >= rajout [
            ask patches with [fixed != TRUE and abs pxcor <= distance_centre and abs pycor <= distance_centre][set landUse 1 set landCover x set pcolor 55 - 10 * x set fixed TRUE]
            set counter counter + rajout
            set distance_centre distance_centre + 1
          ][
            print("final")

            print((word "cb " (objective - counter)))
            ask n-of (objective - counter) patches with [abs pxcor <= distance_centre and abs pycor <= distance_centre and fixed != TRUE][set landUse 1 set landCover x set pcolor 55 - 10 * x set fixed TRUE]
            set counter counter + (objective - counter)
          ]
          print((word "il y a "count patches with [landCover = x] " patches de C"x))
        ]
      ])

    print((word "SNH: " count patches with [landCover = 0]))


  ]

  set totalCiPatches (list count(patches with [landCover = 1]) count(patches with [landCover = 2]) count(patches with [landCover = 3]))
  set proportionCrops (list
    (precision (count(patches with [landCover = 1]) / totalPatchesLandscape * 100) 0)
    (precision (count(patches with [landCover = 2]) / totalPatchesLandscape * 100) 0)
    (precision (count(patches with [landCover = 3]) / totalPatchesLandscape * 100) 0))

  ; pest dynamic variables
  set nbOfCiToInfectPerTick (list)
  set infectionFrequency (list)
  set theoreticalNbFinallyInfectedCrops (list)
  (foreach cropsInfectionrates totalCiPatches  sensibilityframes
    [
      [x y z] ->
      let tempo x * y / 100 ; nombre de patches à infecter à chaque tick
      let lengthframe (item 1 z) - (item 0 z ) + 1
      ifelse tempo = 0
      [
        set nbOfCiToInfectPerTick fput 0 nbOfCiToInfectPerTick
        set infectionFrequency fput round 0 infectionFrequency
        set theoreticalNbFinallyInfectedCrops fput 0 theoreticalNbFinallyInfectedCrops
      ]
      [
       ifelse tempo < 1
       [; on infecte un patch tous les n ticks (défini par infection-pattern-frequency = n)
         set nbOfCiToInfectPerTick fput 1 nbOfCiToInfectPerTick
         set infectionFrequency fput round (1 / tempo) infectionFrequency
         set theoreticalNbFinallyInfectedCrops fput (min(list y (round (lengthframe / round (1 / tempo) )))) theoreticalNbFinallyInfectedCrops
       ]
       [; on infecte n patchs tous les ticks (infection-pattern-frequency = 1)
         set nbOfCiToInfectPerTick fput floor (tempo) nbOfCiToInfectPerTick;; attention ici on arrondi à 0 ça peut amener à un mauvais nb de patches infecté à la fin
         set infectionFrequency fput 1 infectionFrequency
         set theoreticalNbFinallyInfectedCrops fput (min (list y (lengthframe * round( tempo )))) theoreticalNbFinallyInfectedCrops
       ]
      ]

    ] )

  set nbOfCiToInfectPerTick reverse nbOfCiToInfectPerTick
  set infectionFrequency reverse infectionFrequency
  set theoreticalNbFinallyInfectedCrops reverse theoreticalNbFinallyInfectedCrops

  colorPatchesMethod
  computationProximityWithSNHMethod
  findClustersMethod

  ask patches [
   ifelse landCover > 0
   [
     set sensibilityD1 item 0 (item (landCover - 1) sensibilityFrames)
     set sensibilityD2 item 1 (item (landCover - 1) sensibilityFrames)
     set infectionRate item (landCover - 1) cropsInfectionRates
   ]
   [
     set sensibilityD1 0
     set sensibilityD2 0
   ]
  ]

  set totalTicksSimulation (totalYearsSimulation * lengthSeason)



  initializeAdultPredatorsMethod
  set date 0
  set year 1
  set-current-directory (word "/home/tess/Desktop/stage3A/R/traitement_donnees/data/" folderPath)
  initializeFileNamesMethod
end

to go
  if (cropDistribution = "randomly" or cropDistribution = "randomly9" or cropDistribution = "randomly25") and cropsOrder != randomCropsOrder
  [
    show "on s'arrette: cropsOrder != randomCropsOrder"
    stop
  ] ;sert à court-circuiter les combinaisons "random" et les différents crop.order

  ifelse date = 0
  [
    ;;; crops-pests-predators dynamics
    transitionBetweenYearsMethod

    ;;; outputs
    ;writeTickFileMethod
    ;if year = totalyearsSimulation [
      ;ask adultPredators [writeMovementsMethod]]

    show "step1: Year starts"

    ;;; time
    set date date + 1
    tick
  ] ;end if1
  [
    ifelse date <= dateStartForagingForSNH and (count patches with [landUse = 1 and totalCurationCycles > 0] < sum theoreticalNbFinallyInfectedCrops)
    [
      show (word "step2: Foraging for food " date)

      ;;; crops-pests
      ;update-crop-patches-between-start-season-and-date-to-flee
      pestAttackMethod
      switchCropStatusMethod
      updateCropsDuringForagingForFoodMethod
      colorCropPatchesMethod

      ;;; adult-predators dynamic
      mortalityMethod
      foragingForFoodMethod

      ;;; outputs
      ;writeTickFileMethod
      ;if year = totalyearsSimulation [ask adultPredators [writeMovementsMethod]] attention si jamais on veux le fichier adult-predator il y a un appel dans le fichier netlogo adult-predators de netlogo (pour quand ils meurent)



      ;;; time
      tick
      set date date + 1
    ]; end if2
    [
      ifelse date < lengthSeason
      [

        show (word "step3: Foraging for ESN " date)

        ;;; predators of year - 1 die
        ask adultPredators with [birthYear = year - 1][die]


        ;;; crops-pests dynamic
        ifelse date < dateStartForagingForSNH
        [
          ; pestAttackMethod
          switchCropStatusMethod
          updateCropsDuringForagingForFoodMethod
          colorCropPatchesMethod
        ];end if4
        [
          ; pestAttackMethod
          switchCropStatusMethod
          updateCropsDuringForagingForSNHMethod
          colorCropPatchesMethod
        ] ;end else 4

        ;;; adult-predators dynamic
        mortalityMethod
        foragingForSNHMethod

        ;;; outputs
        ;writeTickFileMethod
        ;if year = totalyearsSimulation [ask adultPredators [writeMovementsMethod]]

        ;;; time
        tick
        set date date + 1

      ]; end if 3
      [
        show "step4: Overwintering"

        ;;; kill juveniles and adults in crops because they can not do overwintering
        ask juvenilePredators [ask patch-here [set predatorPresence FALSE] die]
        ask adultPredators with [ [landUse] of patch-here = 1] [ ask patch-here [set predatorPresence FALSE] die ]


        ;;;
        ask patches with [infectionDate = FALSE][set infectionDate 180]
        ;computeRegulationMethod


        ;;; outputs
        ;writeTickFileMethod
        ;writeEndSeasonFileMethod
        if year = TotalYearsSimulation [
          writeSpatialFileMethod
          ;ask adultPredators [writeMovementsMethod]
        ]


        ;;; adult-predators dynamic
        overwinteringMethod ; si on met overwintering ici il sera compté sur l'année n + 1, si on veut en tenir compte en année n, il faut le mettre avant le write

        ;;; indicators
        ;computeRegulationMethod
        ;computeCropLossMethod

        ;;; time : year n -> year n+1
        set year year + 1
        set date 0
      ]; end else3
    ];end else2
  ] ; end else 1

  if ticks >= (totalTicksSimulation + 1) [stop]
end

;to saveVideo [pas]
;
;  vid:start-recorder
;  vid:record-view
;  repeat pas [
;    go
;    vid:record-view
;  ]
;  vid:save-recording "my-movie.mp4"
;
;end

to go_limit [nombre]
  repeat nombre [go]
  ask patches with [infectiondate != false][set pcolor scale-color pink infectiondate 0 50]

  ;ask patches with [timeForCroploss != "NA"and landCover = 3][
    ;ifelse timeForCroploss < 10 [set pcolor 115][
      ;ifelse timeForCroploss < 20 [ set pcolor 125 ][set pcolor 135]]]
  ;show count patches with [pcolor >= 115]
end
@#$#@#$#@
GRAPHICS-WINDOW
540
10
1043
514
-1
-1
15.0
1
10
1
1
1
0
1
1
1
-16
16
-16
16
0
0
1
ticks
30.0

INPUTBOX
11
10
166
70
proportionSNHPatches
7.0
1
0
Number

INPUTBOX
12
141
146
201
proportionCropsPatches
31-31-
1
0
String

TEXTBOX
186
82
336
112
NIL
11
0.0
1

CHOOSER
12
206
150
251
cropDistribution
cropDistribution
"aroundSNH" "randomly" "randomly9" "randomly25"
0

INPUTBOX
153
208
261
268
cropsOrder
2-3-1-
1
0
String

TEXTBOX
273
218
406
250
depuis les SNH vers l'extérieur
11
0.0
1

INPUTBOX
171
10
282
70
targetForAgregation
7.0
1
0
Number

INPUTBOX
286
10
400
70
accuracyThreshold
0.001
1
0
Number

INPUTBOX
403
10
502
70
maxPermutations
500.0
1
0
Number

BUTTON
1053
422
1116
455
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
1053
459
1116
492
NIL
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

TEXTBOX
1057
371
1207
413
C1: jaune\nC2: marron\nC3: orange
11
0.0
1

TEXTBOX
81
362
231
380
NE dynamic
11
0.0
1

INPUTBOX
404
73
517
133
totalYearsSimulation
7.0
1
0
Number

INPUTBOX
102
446
249
506
probaBirthJuvenilePredators
1.0
1
0
Number

INPUTBOX
1121
425
1233
485
folderPath
/variation_foraging/
1
0
String

INPUTBOX
132
383
269
443
dateStartForagingForSNH
150.0
1
0
Number

INPUTBOX
404
137
516
197
lengthSeason
180.0
1
0
Number

INPUTBOX
314
448
376
508
mortality
0.06
1
0
Number

INPUTBOX
277
383
385
443
overwinteringEffect
1.0
1
0
Number

CHOOSER
314
511
425
556
mortalityPattern
mortalityPattern
"constant" "densityDependence"
0

CHOOSER
152
510
310
555
abilityToDetectAttractiveCrops
abilityToDetectAttractiveCrops
9 25 49
0

INPUTBOX
4
382
127
442
foragingMovFrequency
12.0
1
0
Number

CHOOSER
7
509
148
554
foragingPattern
foragingPattern
"Oneof" "RandomB" "Onejump"
0

INPUTBOX
13
271
232
331
framesSensibility
1_50-51_100-101_150-
1
0
String

INPUTBOX
240
271
344
331
infectionRates
40-40-40-
1
0
String

INPUTBOX
7
445
93
505
initNbPredators
100.0
1
0
Number

MONITOR
1056
29
1113
74
NIL
date
17
1
11

MONITOR
1056
131
1116
176
nb crops
length orderCrops
17
1
11

MONITOR
1056
77
1113
122
NIL
year
17
1
11

TEXTBOX
18
331
238
373
ex: 20_33-50_150-60_120-0_150-\n! le 150 c'est récolté donc pas d'attaque possible
11
0.0
1

TEXTBOX
1055
500
1295
598
- dans les notes on prend tjrs l'ex de 4 cultures\n- c'est le pourcentage total du paysage\n- ne pas mettre la dernière, calculé automatiquement pour 100%\nC1-C2-C3-\n20-35-30-
11
0.0
1

TEXTBOX
246
331
461
391
4-10-9-120-\n4 pour 0.4% (c'est divisé par 10 ensuite, permet de s'affranchir des virgules dans le nom de fichier)\n
11
0.0
1

INPUTBOX
12
78
90
138
cropsNumber
3.0
1
0
Number

MONITOR
1412
173
1499
218
NIL
orderCrops
17
1
11

MONITOR
1160
173
1285
218
NIL
cropsInfectionRates
17
1
11

MONITOR
1160
272
1298
317
NIL
nbOfCiToInfectPerTick
17
1
11

MONITOR
1159
223
1252
268
NIL
totalCiPatches
17
1
11

MONITOR
1256
223
1409
268
NIL
sensibilityFrames
17
1
11

MONITOR
1162
320
1281
365
NIL
infectionFrequency
17
1
11

MONITOR
1164
367
1365
412
NIL
theoreticalNbFinallyInfectedCrops
17
1
11

PLOT
1244
10
1444
160
infcrops
NIL
NIL
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"C1" 1.0 0 -1184463 true "" "plot item 0 counterInfectedCrops"
"C2" 1.0 0 -6459832 true "" "plot item 1 counterInfectedCrops"
"C3" 1.0 0 -955883 true "" "plot item 2 counterInfectedCrops"

PLOT
1389
281
1589
431
NE
NIL
NIL
181.0
360.0
0.0
350.0
false
false
"" ""
PENS
"C1" 1.0 0 -10899396 true "" "plot (count adultPredators-on patches with[landCover = 0])"
"pen-1" 1.0 0 -7171555 true "" "plot (count adultPredators-on patches with[landCover = 1])"
"pen-2" 1.0 0 -6459832 true "" "plot (count adultPredators-on patches with[landCover = 2])"
"pen-3" 1.0 0 -955883 true "" "plot (count adultPredators-on patches with[landCover = 3])"
"pen-4" 1.0 0 -7500403 true "" "plot(count adultPredators)"

MONITOR
1289
173
1409
218
NIL
proportionCrops
17
1
11

CHOOSER
98
78
283
123
behaviorSpace
behaviorSpace
"crop_prop_fixed" "crop_prop_variable"
0

INPUTBOX
570
524
683
584
randomCropsOrder
1-2-3-
1
0
String

TEXTBOX
684
542
834
560
pour le behaviorSpace
12
0.0
1

PLOT
1442
520
1642
670
plot 1
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot count(adultpredators)"

BUTTON
927
622
1091
655
NIL
setup2
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

CHOOSER
1092
613
1233
658
config
config
"haut gauche" "centre123" "bande"
1

MONITOR
1100
223
1157
268
SNH nb
count patches with[landCover = 0 ]
17
1
11

INPUTBOX
406
385
515
445
layEgg
7.0
1
0
Number

INPUTBOX
406
447
513
507
juvenileMaturation
14.0
1
0
Number

MONITOR
458
532
558
577
NEProduction
layEgg + juvenileMaturation
17
1
11

PLOT
1516
105
1716
255
NE resource
NIL
NIL
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"C1" 1.0 0 -1184463 true "" "plot count patches with [landCover = 1 and state != 0 ]"
"C2" 1.0 0 -6459832 true "" "plot count patches with [landCover = 2 and state != 0]"
"C3" 1.0 0 -955883 true "" "plot count patches with [landCover = 3 and state != 0]"

TEXTBOX
308
290
347
308
(x10)
12
0.0
1

@#$#@#$#@
## WHAT IS IT?

(a general understanding of what the model is trying to show or explain)

## HOW IT WORKS

(what rules the agents use to create the overall behavior of the model)

## HOW TO USE IT

(how to use the model, including a description of each of the items in the Interface tab)

## THINGS TO NOTICE

(suggested things for the user to notice while running the model)

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

(a reference to the model's URL on the web if it has one, as well as any other necessary credits, citations, and links)
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
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

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

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="foraging freq variation" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="folderPath">
      <value value="&quot;/variation_foraging/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionSNHPatches">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behaviorSpace">
      <value value="&quot;crop_prop_fixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsNumber">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionCropsPatches">
      <value value="&quot;30-30-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropDistribution">
      <value value="&quot;randomly25&quot;"/>
      <value value="&quot;randomly&quot;"/>
      <value value="&quot;randomly9&quot;"/>
      <value value="&quot;aroundSNH&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
      <value value="&quot;3-2-1-&quot;"/>
      <value value="&quot;1-3-2-&quot;"/>
      <value value="&quot;3-1-2-&quot;"/>
      <value value="&quot;2-1-3-&quot;"/>
      <value value="&quot;2-3-1-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectionRates">
      <value value="&quot;40-40-40-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="framesSensibility">
      <value value="&quot;1_50-51_100-101_150-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totalYearsSimulation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lengthSeason">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dateStartForagingForSNH">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPermutations">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accuracyThreshold">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="overwinteringEffect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="probaBirthJuvenilePredators">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initNbPredators">
      <value value="100"/>
    </enumeratedValueSet>
    <steppedValueSet variable="foragingMovFrequency" first="1" step="1" last="12"/>
    <enumeratedValueSet variable="abilityToDetectAttractiveCrops">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingPattern">
      <value value="&quot;Oneof&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortalityPattern">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="randomCropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="infection variation" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="folderPath">
      <value value="&quot;/variation_infection/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionSNHPatches">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behaviorSpace">
      <value value="&quot;crop_prop_fixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsNumber">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionCropsPatches">
      <value value="&quot;30-30-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropDistribution">
      <value value="&quot;aroundSNH&quot;"/>
      <value value="&quot;randomly&quot;"/>
      <value value="&quot;randomly9&quot;"/>
      <value value="&quot;randomly25&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
      <value value="&quot;3-2-1-&quot;"/>
      <value value="&quot;1-3-2-&quot;"/>
      <value value="&quot;3-1-2-&quot;"/>
      <value value="&quot;2-1-3-&quot;"/>
      <value value="&quot;2-3-1-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectionRates">
      <value value="&quot;20-20-20-&quot;"/>
      <value value="&quot;30-30-30-&quot;"/>
      <value value="&quot;40-40-40-&quot;"/>
      <value value="&quot;50-50-50-&quot;"/>
      <value value="&quot;60-60-60-&quot;"/>
      <value value="&quot;70-70-70-&quot;"/>
      <value value="&quot;80-80-80-&quot;"/>
      <value value="&quot;90-90-90-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="framesSensibility">
      <value value="&quot;1_50-51_100-101_150-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totalYearsSimulation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lengthSeason">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dateStartForagingForSNH">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPermutations">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accuracyThreshold">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="overwinteringEffect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="probaBirthJuvenilePredators">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initNbPredators">
      <value value="76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingMovFrequency">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abilityToDetectAttractiveCrops">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingPattern">
      <value value="&quot;Oneof&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortalityPattern">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="randomCropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="overwintering variation" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="folderPath">
      <value value="&quot;/variation_overwintering/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionSNHPatches">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behaviorSpace">
      <value value="&quot;crop_prop_fixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsNumber">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionCropsPatches">
      <value value="&quot;30-30-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropDistribution">
      <value value="&quot;aroundSNH&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
      <value value="&quot;3-2-1-&quot;"/>
      <value value="&quot;1-3-2-&quot;"/>
      <value value="&quot;3-1-2-&quot;"/>
      <value value="&quot;2-1-3-&quot;"/>
      <value value="&quot;2-3-1-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectionRates">
      <value value="&quot;40-40-40-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="framesSensibility">
      <value value="&quot;1_50-51_100-101_150-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totalYearsSimulation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lengthSeason">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dateStartForagingForSNH">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPermutations">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accuracyThreshold">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.06"/>
    </enumeratedValueSet>
    <steppedValueSet variable="overwinteringEffect" first="1" step="1" last="10"/>
    <enumeratedValueSet variable="probaBirthJuvenilePredators">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initNbPredators">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingMovFrequency">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abilityToDetectAttractiveCrops">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingPattern">
      <value value="&quot;Oneof&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortalityPattern">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="randomCropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="mortality variation" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="folderPath">
      <value value="&quot;/variation_mortality/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionSNHPatches">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behaviorSpace">
      <value value="&quot;crop_prop_fixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsNumber">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionCropsPatches">
      <value value="&quot;30-30-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropDistribution">
      <value value="&quot;randomly&quot;"/>
      <value value="&quot;randomly9&quot;"/>
      <value value="&quot;randomly22&quot;"/>
      <value value="&quot;aroundSNH&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
      <value value="&quot;3-2-1-&quot;"/>
      <value value="&quot;1-3-2-&quot;"/>
      <value value="&quot;3-1-2-&quot;"/>
      <value value="&quot;2-1-3-&quot;"/>
      <value value="&quot;2-3-1-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectionRates">
      <value value="&quot;40-40-40-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="framesSensibility">
      <value value="&quot;1_50-51_100-101_150-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totalYearsSimulation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lengthSeason">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dateStartForagingForSNH">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPermutations">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accuracyThreshold">
      <value value="0.001"/>
    </enumeratedValueSet>
    <steppedValueSet variable="mortality" first="0.02" step="0.02" last="0.3"/>
    <enumeratedValueSet variable="overwinteringEffect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="probaBirthJuvenilePredators">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initNbPredators">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingMovFrequency">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abilityToDetectAttractiveCrops">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingPattern">
      <value value="&quot;Oneof&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortalityPattern">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="randomCropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="agregation variation" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="folderPath">
      <value value="&quot;/variation_agregation/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionSNHPatches">
      <value value="7"/>
    </enumeratedValueSet>
    <steppedValueSet variable="targetForAgregation" first="1" step="1" last="9"/>
    <enumeratedValueSet variable="behaviorSpace">
      <value value="&quot;crop_prop_fixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsNumber">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionCropsPatches">
      <value value="&quot;30-30-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropDistribution">
      <value value="&quot;randomly&quot;"/>
      <value value="&quot;randomly9&quot;"/>
      <value value="&quot;randomly22&quot;"/>
      <value value="&quot;aroundSNH&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
      <value value="&quot;3-2-1-&quot;"/>
      <value value="&quot;1-3-2-&quot;"/>
      <value value="&quot;3-1-2-&quot;"/>
      <value value="&quot;2-1-3-&quot;"/>
      <value value="&quot;2-3-1-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectionRates">
      <value value="&quot;40-40-40-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="framesSensibility">
      <value value="&quot;1_50-51_100-101_150-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totalYearsSimulation">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lengthSeason">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dateStartForagingForSNH">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPermutations">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accuracyThreshold">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="overwinteringEffect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="probaBirthJuvenilePredators">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initNbPredators">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingMovFrequency">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abilityToDetectAttractiveCrops">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingPattern">
      <value value="&quot;Oneof&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortalityPattern">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="randomCropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="foraging freq variation random" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="folderPath">
      <value value="&quot;/variation_foraging_aleatoire/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionSNHPatches">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behaviorSpace">
      <value value="&quot;crop_prop_fixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsNumber">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionCropsPatches">
      <value value="&quot;30-30-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropDistribution">
      <value value="&quot;aroundSNH&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
      <value value="&quot;3-2-1-&quot;"/>
      <value value="&quot;1-3-2-&quot;"/>
      <value value="&quot;3-1-2-&quot;"/>
      <value value="&quot;2-1-3-&quot;"/>
      <value value="&quot;2-3-1-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectionRates">
      <value value="&quot;40-40-40-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="framesSensibility">
      <value value="&quot;1_50-51_100-101_150-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totalYearsSimulation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lengthSeason">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dateStartForagingForSNH">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPermutations">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accuracyThreshold">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="overwinteringEffect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="probaBirthJuvenilePredators">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initNbPredators">
      <value value="100"/>
    </enumeratedValueSet>
    <steppedValueSet variable="foragingMovFrequency" first="1" step="1" last="12"/>
    <enumeratedValueSet variable="abilityToDetectAttractiveCrops">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingPattern">
      <value value="&quot;RandomB&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortalityPattern">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="randomCropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="foraging freq variation 1crop et3I" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="folderPath">
      <value value="&quot;/variation_foraging_1crop_et3I/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionSNHPatches">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behaviorSpace">
      <value value="&quot;crop_prop_fixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsNumber">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionCropsPatches">
      <value value="&quot;30-30-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropDistribution">
      <value value="&quot;aroundSNH&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectionRates">
      <value value="&quot;100-&quot;"/>
      <value value="&quot;150-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="framesSensibility">
      <value value="&quot;1_50-51_100-101_150-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totalYearsSimulation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lengthSeason">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dateStartForagingForSNH">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPermutations">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accuracyThreshold">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="overwinteringEffect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="probaBirthJuvenilePredators">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initNbPredators">
      <value value="100"/>
    </enumeratedValueSet>
    <steppedValueSet variable="foragingMovFrequency" first="1" step="1" last="12"/>
    <enumeratedValueSet variable="abilityToDetectAttractiveCrops">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingPattern">
      <value value="&quot;Oneof&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortalityPattern">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="randomCropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="simul1" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="folderPath">
      <value value="&quot;/simul1/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionSNHPatches">
      <value value="2"/>
      <value value="7"/>
      <value value="15"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behaviorSpace">
      <value value="&quot;crop_prop_fixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsNumber">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionCropsPatches">
      <value value="&quot;30-30-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropDistribution">
      <value value="&quot;randomly22&quot;"/>
      <value value="&quot;randomly&quot;"/>
      <value value="&quot;randomly9&quot;"/>
      <value value="&quot;aroundSNH&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
      <value value="&quot;3-2-1-&quot;"/>
      <value value="&quot;1-3-2-&quot;"/>
      <value value="&quot;3-1-2-&quot;"/>
      <value value="&quot;2-1-3-&quot;"/>
      <value value="&quot;2-3-1-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="randomCropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectionRates">
      <value value="&quot;40-40-40-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="framesSensibility">
      <value value="&quot;1_50-51_100-101_150-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totalYearsSimulation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lengthSeason">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dateStartForagingForSNH">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPermutations">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accuracyThreshold">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="overwinteringEffect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="probaBirthJuvenilePredators">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initNbPredators">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingMovFrequency">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abilityToDetectAttractiveCrops">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingPattern">
      <value value="&quot;Oneof&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortalityPattern">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="dispersal variation" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="folderPath">
      <value value="&quot;/variation_dispersal/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionSNHPatches">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behaviorSpace">
      <value value="&quot;crop_prop_fixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsNumber">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionCropsPatches">
      <value value="&quot;30-30-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropDistribution">
      <value value="&quot;aroundSNH&quot;"/>
      <value value="&quot;randomly&quot;"/>
      <value value="&quot;randomly9&quot;"/>
      <value value="&quot;randomly22&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
      <value value="&quot;3-2-1-&quot;"/>
      <value value="&quot;1-3-2-&quot;"/>
      <value value="&quot;3-1-2-&quot;"/>
      <value value="&quot;2-1-3-&quot;"/>
      <value value="&quot;2-3-1-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectionRates">
      <value value="&quot;40-40-40-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="framesSensibility">
      <value value="&quot;1_50-51_100-101_150-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totalYearsSimulation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lengthSeason">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dateStartForagingForSNH">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPermutations">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accuracyThreshold">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="overwinteringEffect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="probaBirthJuvenilePredators">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initNbPredators">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingMovFrequency">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abilityToDetectAttractiveCrops">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingPattern">
      <value value="&quot;RandomB&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortalityPattern">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="randomCropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="infection variation faible" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="folderPath">
      <value value="&quot;/variation_infection_faible/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionSNHPatches">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behaviorSpace">
      <value value="&quot;crop_prop_fixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsNumber">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionCropsPatches">
      <value value="&quot;30-30-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropDistribution">
      <value value="&quot;aroundSNH&quot;"/>
      <value value="&quot;randomly&quot;"/>
      <value value="&quot;randomly9&quot;"/>
      <value value="&quot;randomly25&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
      <value value="&quot;3-2-1-&quot;"/>
      <value value="&quot;1-3-2-&quot;"/>
      <value value="&quot;3-1-2-&quot;"/>
      <value value="&quot;2-1-3-&quot;"/>
      <value value="&quot;2-3-1-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectionRates">
      <value value="&quot;1-1-1-&quot;"/>
      <value value="&quot;3-3-3-&quot;"/>
      <value value="&quot;5-5-5-&quot;"/>
      <value value="&quot;8-8-8-&quot;"/>
      <value value="&quot;1O-10-10-&quot;"/>
      <value value="&quot;15-15-15-&quot;"/>
      <value value="&quot;20-20-20-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totalYearsSimulation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lengthSeason">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dateStartForagingForSNH">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPermutations">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accuracyThreshold">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="overwinteringEffect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="probaBirthJuvenilePredators">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initNbPredators">
      <value value="76"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingMovFrequency">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abilityToDetectAttractiveCrops">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingPattern">
      <value value="&quot;Oneof&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortalityPattern">
      <value value="&quot;constant&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="randomCropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="simul3_test" repetitions="10" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <enumeratedValueSet variable="folderPath">
      <value value="&quot;/tests/&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionSNHPatches">
      <value value="2"/>
      <value value="7"/>
      <value value="15"/>
      <value value="30"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="targetForAgregation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="behaviorSpace">
      <value value="&quot;crop_prop_fixed&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsNumber">
      <value value="3"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="proportionCropsPatches">
      <value value="&quot;30-30-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropDistribution">
      <value value="&quot;randomly25&quot;"/>
      <value value="&quot;randomly&quot;"/>
      <value value="&quot;randomly9&quot;"/>
      <value value="&quot;aroundSNH&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
      <value value="&quot;3-2-1-&quot;"/>
      <value value="&quot;1-3-2-&quot;"/>
      <value value="&quot;3-1-2-&quot;"/>
      <value value="&quot;2-1-3-&quot;"/>
      <value value="&quot;2-3-1-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="randomCropsOrder">
      <value value="&quot;1-2-3-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="infectionRates">
      <value value="&quot;10-10-10-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="framesSensibility">
      <value value="&quot;1_50-51_100-101_150-&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="totalYearsSimulation">
      <value value="7"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="lengthSeason">
      <value value="180"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="dateStartForagingForSNH">
      <value value="150"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="maxPermutations">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="accuracyThreshold">
      <value value="0.001"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortality">
      <value value="0.06"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="overwinteringEffect">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="probaBirthJuvenilePredators">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initNbPredators">
      <value value="100"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingMovFrequency">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="abilityToDetectAttractiveCrops">
      <value value="9"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="foragingPattern">
      <value value="&quot;Oneof&quot;"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="mortalityPattern">
      <value value="&quot;constant&quot;"/>
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
