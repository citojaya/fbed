#!MC 900
$!VarSet |MFBD| = 'C:\Program Files\TEC90'

### MACRO DEFINITION ###
$!VARSET |ROOT_FILE| = "E:\Documents\fbed\wall"
$!VARSET |INPUT_FILE| = "|ROOT_FILE|.PLT"
$!VARSET |OUTPUT_FILE| = "|ROOT_FILE|_ooo.avi"

###### total zone: #######
$!VARSET |CHAMBER_ZONE| = 1
$!VARSET |START_ZONE| = 2
$!VARSET |END_ZONE| = 60
$!VARSET |ZONE_INCREMENT| = 1

$!VARSET |START_TIME| = 1.5
$!VARSET |TIME_INTERVAL| = 0.001

$!VARSET |AVI_SPEED| = 10


### READ INPUT DATA FILE ###
$!READDATASET  "|INPUT_FILE|"
  READDATAOPTION = NEW
  RESETSTYLE = YES
  INCLUDETEXT = NO
  INCLUDEGEOM = NO
  INCLUDECUSTOMLABELS = NO
  VARLOADMODE = BYNAME
  VARNAMELIST = '"V1" "V2" "V3" "V4" "V5" "V6" "V7" "V8"'
  ZONELIST = [|CHAMBER_ZONE|,|START_ZONE|-|END_ZONE|]

### OUTPUT FILE ###
$!EXPORTSETUP
  EXPORTFORMAT = AVI
  ANIMATIONSPEED = |AVI_SPEED|
  IMAGEWIDTH = 880
  EXPORTFNAME = "|OUTPUT_FILE|"

### TEXT OPTION ###
$!ATTACHTEXT
  XYPOS
    {
    X = 40
    Y = 10
    }
  COLOR = BLACK
  TEXTSHAPE
    {
    FONT = HELVBOLD
    HEIGHT = 24
    }
TEXT = ' Fludized Bed '
$!ATTACHTEXT
  XYPOS
    {
    X = 79
    Y = 32
    }
  COLOR = BLACK
  TEXTSHAPE
    {
    FONT = HELVBOLD
    HEIGHT = 14
    }
TEXT = ' Velocity (m/s) '

### COLOR OPTION ###
$!COLORMAP
  CONTOURCOLORMAP = SMRAINBOW
$!COLORMAPCONTROL RESETTOFACTORY

### Frame Number 1 ###
$!FRAMELAYOUT
  SHOWBORDER = NO
  BACKGROUNDCOLOR = WHITE
  XYPOS
  {
    X = -0.013293
    Y = 0.029154
  }
  WIDTH = 10
  HEIGHT = 8
$!FRAMEMODE  = THREED
$!FRAMENAME  = 'Frame 001'

### CONTOUR AND LENGEND OPTION ###
$!GLOBALCONTOUR
  VAR = 8
  LABELS
  {
    AUTOLEVELSKIP = 2
  }
  LEGEND
  {
    SHOW = YES
    SHOWHEADER = NO
    TEXTCOLOR = BLACK
    XYPOS
    {
      X = 90
      Y = 70
    }
    OVERLAYBARGRID = NO
    NUMFORMAT
    {
      FORMATTING = FIXEDFLOAT
      PRECISION = 1
    }
  }
  COLORCUTOFF
  {
    RANGEMIN = 0
    RANGEMAX = 30
  }
  COLORMAPFILTER
  {
    CONTINUOUSCOLOR
    {
      CMAX = 30
    }
  }
  
$!CONTOURLEVELS NEW
  RAWDATA
20
0
0.15
0.3
0.45
0.6
0.75
0.9
1.05
1.20
1.35
1.50
2.65
2.80
2.95
4.10
4.25
4.35
4.50
4.65
20


### GLOBAL SCATTER OPTION ###
$!GLOBALSCATTER
  VAR = 7
  RELATIVESIZE = 1
  
  
### AXIS OPTION ###
$!THREEDAXIS
  XVAR = 1
  YVAR = 2
  ZVAR = 3
  FRAMEAXIS
  {
  	SHOW = YES
  	COLOR = BLACK
  	XYPOS
    {
      X = 8.9637
      Y = 20.46
    }
  }
### $!VIEW FIT ###
$!THREEDAXIS
  AXISMODE = XYZDEPENDENT
  XYDEPXTOYRATIO = 1
  DEPXTOYRATIO = 1
  DEPXTOZRATIO = 1
  XEDGE = 3
  YEDGE = 3
  XDETAIL
    {
    SHOWAXIS = YES
    RANGEMIN = -18.6
    RANGEMAX = 12.6
    GRSPACING = 10
    }
  YDETAIL
    {
    SHOWAXIS = YES
    RANGEMIN = -16.6
    RANGEMAX = 16.6
    GRSPACING = 10
    }
  ZDETAIL
    {
    SHOWAXIS = YES
    RANGEMIN = -4.16875
    RANGEMAX = 70.79374999999
    GRSPACING = 5
    }
$!GLOBALISOSURFACE
  ISOVALUE1 = 5.61518001556
  ISOVALUE2 = 11.2303600311
  ISOVALUE3 = 16.8455400467
  
$!GLOBALTHREED
  AXISSCALEFACT
  {
    X = 1
    Y = 1
    Z = 1
  }
  ROTATEORIGIN
  {
    X = 0
    Y = 0
    Z = 1.8125
  }
  LINELIFTFRACTION = 0.5
  SYMBOLLIFTFRACTION = 0.6
  VECTORLIFTFRACTION = 0.7
  
 $!THREEDVIEW
  PSIANGLE = -90  #Rotate about Y axis
  THETAANGLE = -90
  ALPHAANGLE = 90
  VIEWERPOSITION
    {
    X = -10.84143494188
    Y = 80.52950932747
    Z = 40.791412956
    }
  VIEWWIDTH = 360

#$!THREEDVIEW
# PSIANGLE = 22.349
 # THETAANGLE = -118.869
  #ALPHAANGLE = 116.767
  #VIEWERPOSITION
   # {
    #X = 95.535136064
    #Y = 56.9442313489
    #Z = 285.531513957
    #}
  #VIEWWIDTH = 200
  
$!FIELDLAYERS
  SHOWMESH = NO
  SHOWSCATTER = YES
  SHOWSHADE = YES
  SHOWBOUNDARY = NO

# ANIMATION START HERE

$!VARSET |ZONE1| = |START_ZONE|
$!VARSET |ZONE2| = |ZONE1|
$!VARSET |ZONE2| += 1

$!VARSET |SIMTIME| = |START_TIME|

$!ACTIVEFIELDZONES = [|CHAMBER_ZONE|]
#$!ACTIVEFIELDZONES = [|ZONE1|-|ZONE2|]
$!FIELD [|CHAMBER_ZONE|]
	MESH
	{
	  SHOW = NO
	}
	SCATTER
	{
		SHOW = NO
	}
	SHADE
	{
		SHOW = NO
		SHADETYPE = COLOREDPANELED
		COLOR = BLUE
	}

$!VARSET |ZONE1| = |START_ZONE|
$!VARSET |ZONE2| = |ZONE1|
$!VARSET |ZONE2| += 1

### STARTING LOOP #######

### LOOP_NUM = (END_ZONE - START_ZONE + 1)/2
$!VARSET |LOOP_NUM| = |END_ZONE|
$!VARSET |LOOP_NUM| -= |START_ZONE|
$!VARSET |LOOP_NUM| += 1
$!VARSET |LOOP_NUM| /= |ZONE_INCREMENT|

$!LOOP |LOOP_NUM|


$!ACTIVEFIELDZONES += [|ZONE1|]
$!FIELD [|ZONE1|]
	MESH
  {
    SHOW = NO
  }

	SCATTER
  {
    SHOW = NO
  }
	SHADE
	{	SHOW =YES
		SHADETYPE = COLOREDPANELED
		COLOR = CUSTOM7
	}

$!ACTIVEFIELDZONES += [|ZONE2|]
$!FIELD [|ZONE2|]
	MESH
  {
    SHOW = NO
  }
	SCATTER
	{
		SHOW = YES
		SYMBOLSHAPE
		{
			ISASCII = NO
			GEOMSHAPE = CIRCLE
		}
		COLOR = BLACK
		ISFILLED = YES
		FILLCOLOR = MULTI
		SIZEBYVARIABLE = YES
	}
	SHADE
	{
	  SHOW = NO
	}

$!ATTACHTEXT
  XYPOS
    {
    X = 40
    Y = 93
    }
  COLOR = BLACK
  BOX
    {
    BOXTYPE = FILLED
    COLOR = WHITE
    FILLCOLOR = WHITE
    }
  TEXTSHAPE
    {
    FONT = HELVBOLD
    HEIGHT = 28
    }

  #TEXT = 'Cross section at middle part of the mill'
  TEXT = 't = |SIMTIME%6.3f| sec'

#$!BLANKING IJK{ZONE = 1}
#$!BLANKING IJK{INCLUDE = YES}
#$!BLANKING IJK{IJKBLANKMODE = INTERIOR}
#$!BLANKING IJK{IJKBLANKMODE = EXTERIOR}
#$!BLANKING IJK{JMAXFRACT = 2}
#$!BLANKING IJK{KMINFRACT = 2}
#$!BLANKING IJK{KMAXFRACT = 5}

#$!BLANKING VALUE{INCLUDE = YES}
#$!BLANKING VALUE{CONSTRAINT 1 {INCLUDE = YES}}
#$!BLANKING VALUE{CONSTRAINT 1 {VARA = 1}}
#$!BLANKING VALUE{CONSTRAINT 1 {RELOP = GREATERTHANOREQUAL}}
#$!BLANKING VALUE{CONSTRAINT 1 {VALUECUTOFF = 45}}

#$!BLANKING VALUE{INCLUDE = YES}
#$!BLANKING VALUE{CONSTRAINT 1 {INCLUDE = YES}}
#$!BLANKING VALUE{CONSTRAINT 1 {VARA = 3}}
#$!BLANKING VALUE{CONSTRAINT 1 {RELOP = GREATERTHANOREQUAL}}
#$!BLANKING VALUE{CONSTRAINT 1 {VALUECUTOFF = 54}}

$!REDRAWALL
$!IF |LOOP| == 1
	$!EXPORTSTART
$!ENDIF
$!IF |LOOP| != 1
	$!EXPORTNEXTFRAME
$!ENDIF

$!ACTIVEFIELDZONES -= [|ZONE1|]
$!ACTIVEFIELDZONES -= [|ZONE2|]

$!VARSET |ZONE1| += |ZONE_INCREMENT|
$!VARSET |ZONE2| += |ZONE_INCREMENT|

$!VARSET |SIMTIME| += |TIME_INTERVAL|
$!ENDLOOP


$!EXPORTFINISH

$!RemoveVar |MFBD|