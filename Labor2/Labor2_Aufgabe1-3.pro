CoDeSys+'   �                   @        @   2.3.9.31    @/    @                             Js�R +    @          �=              ��MR        �   @   q   C:\TwinCAT\PLC\LIB\STANDARD.LIB @                                                                                          CONCAT               STR1               ��              STR2               ��                 CONCAT                                         ��66  �   ����           CTD           M             ��           Variable for CD Edge Detection      CD            ��           Count Down on rising edge    LOAD            ��           Load Start Value    PV           ��           Start Value       Q            ��           Counter reached 0    CV           ��           Current Counter Value             ��66  �   ����           CTU           M             ��            Variable for CU Edge Detection       CU            ��       
    Count Up    RESET            ��           Reset Counter to 0    PV           ��           Counter Limit       Q            ��           Counter reached the Limit    CV           ��           Current Counter Value             ��66  �   ����           CTUD           MU             ��            Variable for CU Edge Detection    MD             ��            Variable for CD Edge Detection       CU            ��	       
    Count Up    CD            ��
           Count Down    RESET            ��           Reset Counter to Null    LOAD            ��           Load Start Value    PV           ��           Start Value / Counter Limit       QU            ��           Counter reached Limit    QD            ��           Counter reached Null    CV           ��           Current Counter Value             ��66  �   ����           DELETE               STR               ��              LEN           ��              POS           ��                 DELETE                                         ��66  �   ����           F_TRIG           M             ��
                 CLK            ��           Signal to detect       Q            ��           Edge detected             ��66  �   ����           FIND               STR1               ��              STR2               ��                 FIND                                     ��66  �   ����           INSERT               STR1               ��              STR2               ��              POS           ��                 INSERT                                         ��66  �   ����           LEFT               STR               ��              SIZE           ��                 LEFT                                         ��66  �   ����           LEN               STR               ��                 LEN                                     ��66  �   ����           MID               STR               ��              LEN           ��              POS           ��                 MID                                         ��66  �   ����           R_TRIG           M             ��
                 CLK            ��           Signal to detect       Q            ��           Edge detected             ��66  �   ����           REPLACE               STR1               ��              STR2               ��              L           ��              P           ��                 REPLACE                                         ��66  �   ����           RIGHT               STR               ��              SIZE           ��                 RIGHT                                         ��66  �   ����           RS               SET            ��              RESET1            ��                 Q1            ��
                       ��66  �   ����           SEMA           X             ��                 CLAIM            ��	              RELEASE            ��
                 BUSY            ��                       ��66  �   ����           SR               SET1            ��              RESET            ��                 Q1            ��	                       ��66  �   ����           TOF           M             ��           internal variable 	   StartTime            ��           internal variable       IN            ��       ?    starts timer with falling edge, resets timer with rising edge    PT           ��           time to pass, before Q is set       Q            ��	       2    is FALSE, PT seconds after IN had a falling edge    ET           ��
           elapsed time             ��66  �   ����           TON           M             ��           internal variable 	   StartTime            ��           internal variable       IN            ��       ?    starts timer with rising edge, resets timer with falling edge    PT           ��           time to pass, before Q is set       Q            ��	       0    is TRUE, PT seconds after IN had a rising edge    ET           ��
           elapsed time             ��66  �   ����           TP        	   StartTime            ��           internal variable       IN            ��       !    Trigger for Start of the Signal    PT           ��       '    The length of the High-Signal in 10ms       Q            ��	           The pulse    ET           ��
       &    The current phase of the High-Signal             ��66  �   ����    R    @                                                                                          MAIN     
      yrect                             ysin                             tim                            MyArray   	  c                                         iter                            Mean                             Mean_2                          	   SampleNow               	              ActualN              
              Sum_1                                              ��MR  @    ����           MEANERV               Sum            &               N           &                  MEANERV                                      ��MR  @    ����           MEANV           Sum             %               i            %                  YARRAY   	  c                        %                  MEANV                                      ��MR  @    ����           RECTGEN           Output           0    #                  Ymin            #               Ymax            #               Yoff            #               w            #               t            #                  RECTGEN                                      ��MR  @    ����           SINGEN               A            !               Yoff            !               w            !               t            !                  SINGEN                                      ��MR  @    ����            
 �   &   !   #      %       ( �      K   �     K   �     K   �     K   �                 	         +     ��localhost    �:"��  9    &   �9�9�9<� �� 7Dw             H7Dw  O       ɯ            ɯ      ě� X`����0p6� D    F  �� �� �r� ����    x�h4� �}#         ��3�(�       ɯ           ɯ  � ě� � �� ě� p`������ �F�     ,   ,                                                        K         @   ��MR�  /*BECKCONFI3*/
        !�� @   @   �   �     3               
   Standard            	��MR                        VAR_GLOBAL
END_VAR
                                                                                  "   , � � v�             Standard
         MAIN����               ��MR                 $����                                           Standard ��MR	��MR                                   	��MR                        VAR_CONFIG
END_VAR
                                                                                   '                  �              Global_Variables ��MR	��MR       �               VAR_GLOBAL
END_VAR
                                                                                               '           "        ����           TwinCAT_Configuration ��MR	��MR"                     K   (* Generated automatically by TwinCAT - (read only) *)
VAR_CONFIG
END_VAR                                                                                               '           	     ) 1 9              Variable_Configuration ��MR	��MR	       �               VAR_CONFIG
END_VAR
                                                                                                 �   |0|0 @|    @Z   MS Sans Serif @       HH':'mm':'ss @      dd'-'MM'-'yyyy   dd'-'MM'-'yyyy HH':'mm':'ss�����                               4     �   ���  �3 ���   � ���     
    @��  ���     @      DEFAULT             System      �   |0|0 @|    @Z   MS Sans Serif @       HH':'mm':'ss @      dd'-'MM'-'yyyy   dd'-'MM'-'yyyy HH':'mm':'ss�����                      )   HH':'mm':'ss @                             dd'-'MM'-'yyyy @       '         , 4 4 ��           MAIN ӵMR	��MR                      �   PROGRAM MAIN
VAR
	yrect, ysin : REAL;
	tim: INT;
	MyArray: ARRAY[0..99] OF REAL;
	iter: INT;
	Mean: REAL;
	Mean_2: REAL;
	SampleNow: REAL;
	ActualN: INT;
	Sum_1:REAL;
END_VARj  yrect := RECTGEN(Ymin:=-1,Ymax:=1,Yoff:=2,w:=10,t:=tim);
ysin := SINGEN(A:=1,Yoff:=0, w:=10, t:=tim);
tim := (tim+1);

iter := iter+1;
IF ( iter = 100) THEN
iter:= 0;
END_IF;

MyArray[iter]:=yrect;
SampleNow := yrect;

ActualN:= ActualN + 1;
Sum_1 := Sum_1+SampleNow;

Mean := MEANV(YARRAY:= MyArray);
Mean_2 := MEANERV(Sum:= Sum_1, N:= ActualN);               &   , � � \T           MEANERV ��MR	��MR                      P   FUNCTION MEANERV : REAL
VAR_INPUT
	Sum: REAL;
	N: INT;
END_VAR
VAR
END_VAR   MEANERV:= Sum/N;               %   , h h (            MEANV ��MR	��MR                      n   FUNCTION MEANV : REAL
VAR_INPUT
	YARRAY: ARRAY[0..99] OF REAL;
END_VAR
VAR
	Sum: REAL;
	i: INT;
END_VARX   FOR i:= 0 TO 99 BY 1 DO
	Sum := Sum +  YARRAY[i];
END_FOR;
         MEANV := Sum/100;               #   , h h (            RECTGEN ��MR	��MR                      n   FUNCTION RECTGEN : REAL
VAR_INPUT
	Ymin, Ymax, Yoff, w, t: REAL;
END_VAR
VAR
	Output: REAL := 0;
END_VARi   IF ((1*SIN(w*t/100)) < 0) THEN
	Output:= Ymin;
ELSE
	Output:= Ymax;
END_IF;

RECTGEN:= Output+Yoff;               !   , N N            SINGEN ��MR	��MR                      O   FUNCTION SINGEN :REAL
VAR_INPUT
	A, Yoff, w, t:  REAL;
END_VAR
VAR
END_VAR&   SINGEN := A * SIN( w * t/100 ) + Yoff;                 ����                   "   STANDARD.LIB 5.6.98 11:03:02 @V�w5      CONCAT @                	   CTD @        	   CTU @        
   CTUD @           DELETE @           F_TRIG @        
   FIND @           INSERT @        
   LEFT @        	   LEN @        	   MID @           R_TRIG @           REPLACE @           RIGHT @           RS @        
   SEMA @           SR @        	   TOF @        	   TON @           TP @              Global Variables 0 @                          pz��               2                ����������������  
             ����                  ����  wn	Sdece                      POUs                MAIN                      MEANERV  &                   MEANV  %                   RECTGEN  #                   SINGEN  !   ����          
   Data types  ����             Visualizations  ����               Global Variables                 Global_Variables                     TwinCAT_Configuration  "                   Variable_Configuration  	   ����                                                             ��MR                         	   localhost            P      	   localhost            P      	   localhost            P            ���