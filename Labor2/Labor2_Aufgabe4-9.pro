CoDeSys++   �                   @        @   2.3.9.31    @/    @                             4s�R +    @          �=              ��MR        �   @   q   C:\TwinCAT\PLC\LIB\STANDARD.LIB @                                                                                          CONCAT               STR1               ��              STR2               ��                 CONCAT                                         ��66  �   ����           CTD           M             ��           Variable for CD Edge Detection      CD            ��           Count Down on rising edge    LOAD            ��           Load Start Value    PV           ��           Start Value       Q            ��           Counter reached 0    CV           ��           Current Counter Value             ��66  �   ����           CTU           M             ��            Variable for CU Edge Detection       CU            ��       
    Count Up    RESET            ��           Reset Counter to 0    PV           ��           Counter Limit       Q            ��           Counter reached the Limit    CV           ��           Current Counter Value             ��66  �   ����           CTUD           MU             ��            Variable for CU Edge Detection    MD             ��            Variable for CD Edge Detection       CU            ��	       
    Count Up    CD            ��
           Count Down    RESET            ��           Reset Counter to Null    LOAD            ��           Load Start Value    PV           ��           Start Value / Counter Limit       QU            ��           Counter reached Limit    QD            ��           Counter reached Null    CV           ��           Current Counter Value             ��66  �   ����           DELETE               STR               ��              LEN           ��              POS           ��                 DELETE                                         ��66  �   ����           F_TRIG           M             ��
                 CLK            ��           Signal to detect       Q            ��           Edge detected             ��66  �   ����           FIND               STR1               ��              STR2               ��                 FIND                                     ��66  �   ����           INSERT               STR1               ��              STR2               ��              POS           ��                 INSERT                                         ��66  �   ����           LEFT               STR               ��              SIZE           ��                 LEFT                                         ��66  �   ����           LEN               STR               ��                 LEN                                     ��66  �   ����           MID               STR               ��              LEN           ��              POS           ��                 MID                                         ��66  �   ����           R_TRIG           M             ��
                 CLK            ��           Signal to detect       Q            ��           Edge detected             ��66  �   ����           REPLACE               STR1               ��              STR2               ��              L           ��              P           ��                 REPLACE                                         ��66  �   ����           RIGHT               STR               ��              SIZE           ��                 RIGHT                                         ��66  �   ����           RS               SET            ��              RESET1            ��                 Q1            ��
                       ��66  �   ����           SEMA           X             ��                 CLAIM            ��	              RELEASE            ��
                 BUSY            ��                       ��66  �   ����           SR               SET1            ��              RESET            ��                 Q1            ��	                       ��66  �   ����           TOF           M             ��           internal variable 	   StartTime            ��           internal variable       IN            ��       ?    starts timer with falling edge, resets timer with rising edge    PT           ��           time to pass, before Q is set       Q            ��	       2    is FALSE, PT seconds after IN had a falling edge    ET           ��
           elapsed time             ��66  �   ����           TON           M             ��           internal variable 	   StartTime            ��           internal variable       IN            ��       ?    starts timer with rising edge, resets timer with falling edge    PT           ��           time to pass, before Q is set       Q            ��	       0    is TRUE, PT seconds after IN had a rising edge    ET           ��
           elapsed time             ��66  �   ����           TP        	   StartTime            ��           internal variable       IN            ��       !    Trigger for Start of the Signal    PT           ��       '    The length of the High-Signal in 10ms       Q            ��	           The pulse    ET           ��
       &    The current phase of the High-Signal             ��66  �   ����    R    @                                                                                          ERROR                             G�VR  @    ����           INIT                             G�VR  @    ����           MAIN           Modus               Betriebs_Modus                    Bit_0                             Bit_1                             Bit_2                             Bit_3                             Bit_4                             Bit_5               	              Bit_6               
              Bit_7                             Byte_0                         	   Byte_temp                            Word_0                            spannung1INT                           spannung1Real                             PID1                     PID_Parameter                    PID2                     PID_Parameter                    Start                             Stop                             Schwelle                             Speed                             TEST                            Gaggi                                             2�VR  @   ����           RUN                             G�VR  @    ����           SPANNUNGINT_TO_REAL               spannungINT                             SpannungINT_to_REAL                                      G�VR  @    ����           STOPPED                             G�VR  @    ����            
 �      *      '   (   )   #   "      !       ( �      K   �     K   �     K   �     K   �                 �         +     ��localhost    �:"��  9    &   �9�9�9<� �� 7Dw             H7Dw  O       ɯ            ɯ      ě� X`����0p6� D    F  �� �� �r� ����    x�h4� �}#         ��3�(�       ɯ           ɯ  � ě� � �� ě� p`������ �F�     ,   ,                                                        K         @   G�VR�  /*BECKCONFI3*/
        !�� @   @   �   �     3               
   Standard            	2�VR                        VAR_GLOBAL
END_VAR
                                                                                  "   , � � v�             Standard
         MAIN����               2�VR                 $����                                           Standard ��MR	��MR                                   	2�VR                        VAR_CONFIG
END_VAR
                                                                                   '                  �              Global_Variables G�VR	G�VR       �               VAR_GLOBAL
END_VAR
                                                                                               '           "   ,   �           TwinCAT_Configuration G�VR	2�VR"                     �   (* Generated automatically by TwinCAT - (read only) *)
VAR_CONFIG
	MAIN.spannung1INT AT %IB0 : INT;	(*  ~ {LinkedWith:TIID^Ger�t 1 (EtherCAT)^Klemme 1 (EK1200)^Klemme 2 (EL3004)^AI Standard Channel 1^Value} *)
END_VAR                                                                                               '           	     ) 1 9              Variable_Configuration G�VR	G�VR	       �               VAR_CONFIG
END_VAR
                                                                                                 �   |0|0 @|    @Z   MS Sans Serif @       HH':'mm':'ss @      dd'-'MM'-'yyyy   dd'-'MM'-'yyyy HH':'mm':'ss�����                               4     �   ���  �3 ���   � ���     
    @��  ���     @      DEFAULT             System      �   |0|0 @|    @Z   MS Sans Serif @       HH':'mm':'ss @      dd'-'MM'-'yyyy   dd'-'MM'-'yyyy HH':'mm':'ss�����                      )   HH':'mm':'ss @                             dd'-'MM'-'yyyy @       '   *   , � � `x           Betriebs_Modus G�VR	G�VR          �         V   TYPE Betriebs_Modus :
(
	Initial := 0,
	Rn := 1,
	Spt:= 2,
	Err :=3
);
END_TYPE             !   , K K H�           PID_Parameter G�VR	G�VR      	Wd_WO;        �   TYPE PID_Parameter :
STRUCT
	Kp:REAL;
	Ti:REAL;
	Td:REAL;
	N:UINT;
	ARW:BOOL;
	umin:REAL;
	umax:REAL;
END_STRUCT
END_TYPE              )   , � � G_           ERROR G�VR	G�VR      �                 PROGRAM ERROR
VAR
END_VAR   ;                  , K K ��           INIT G�VR	G�VR      �wwwwww           PROGRAM INIT
VAR
END_VAR   ;                   , 4 4 ��           MAIN 2�VR	2�VR                      �  PROGRAM MAIN
VAR
	Modus : Betriebs_Modus;
	Bit_0 : BOOL;
	Bit_1:BOOL;
	Bit_2:BOOL;
	Bit_3:BOOL;
	Bit_4:BOOL;
	Bit_5:BOOL;
	Bit_6:BOOL;
	Bit_7:BOOL;
	Byte_0:BYTE;
	Byte_temp:BYTE;
	Word_0:WORD;
	spannung1INT AT %I*:INT;
	spannung1Real:REAL;
	PID1:PID_Parameter;
	PID2:PID_Parameter;
	Start:BOOL;
	Stop:BOOL;
	Schwelle:REAL;
	Speed:REAL;
END_VAR
VAR PERSISTENT
	TEST: REAL;
END_VAR
VAR RETAIN
	Gaggi: REAL;
END_VARp  CASE Modus OF
Initial:
	INIT();
Rn:
	RUN();
Spt:
	STOPPED();
Err:
	ERROR();
ELSE
	STOPPED();
END_CASE;

Byte_0:=BOOL_TO_BYTE(Bit_0);
Byte_0:=SHL(BOOL_TO_BYTE(Bit_1),1);
Byte_0:=SHL(BOOL_TO_BYTE(Bit_2),2);
Byte_0:=SHL(BOOL_TO_BYTE(Bit_3),3);
Byte_0:=SHL(BOOL_TO_BYTE(Bit_4),4);
Byte_0:=SHL(BOOL_TO_BYTE(Bit_5),5);
Byte_0:=SHL(BOOL_TO_BYTE(Bit_6),6);
Byte_0:=SHL(BOOL_TO_BYTE(Bit_7),7);

Word_0:=Byte_0;
Byte_temp:=16#F;
Word_0:=Word_0 OR SHL(Byte_temp,8);
spannung1Real :=SpannungINT_to_REAL(spannung1INT);


PID1.ARW:=0;
PID1.Kp:=10;
PID1.Td:=0.01;
PID1.Ti:=0.001;

Gaggi:=10;

TEST:=10;
               '   , } } -           RUN G�VR	G�VR      �                  PROGRAM RUN
VAR
END_VAR   ;                  , 2 2 /�           SpannungINT_to_REAL ��VR	G�VR                    X   FUNCTION SpannungINT_to_REAL : REAL
VAR_INPUT
	spannungINT:INT;
END_VAR
VAR
END_VAR9   SpannungINT_to_REAL := INT_TO_REAL(spannungINT)/32768*10;               (   , � � .F           STOPPED G�VR	G�VR      �                  PROGRAM STOPPED
VAR
END_VAR   ;                #   , d d a�           myVisual G�VR
    @    G�VR   d                                                                                                        
    @        2 < � e i P     @                    Start @���     ���             @         ���        
   MAIN.Start                 @       �                                                                                                     
    @        2 x � � i �     @                    Stop @���     ���             @        ���        	   MAIN.Stop                 @       �                                                                                                     
    @        2 � � � � 
   MAIN.Speed   0   10.9                                                                                                            
    @        2 � �    MAIN.Schwelle   -10   15.4                �   ��   �   ��   � � � ���     �   ��   �   ��   � � � ���                  ����                   "   STANDARD.LIB 5.6.98 11:03:02 @V�w5      CONCAT @                	   CTD @        	   CTU @        
   CTUD @           DELETE @           F_TRIG @        
   FIND @           INSERT @        
   LEFT @        	   LEN @        	   MID @           R_TRIG @           REPLACE @           RIGHT @           RS @        
   SEMA @           SR @        	   TOF @        	   TON @           TP @              Global Variables 0 @                          pz��               2                ����������������  
             ����                  ����  wn	Sdece                      POUs                 ERROR  )                   INIT                    MAIN                      RUN  '                   SpannungINT_to_REAL                     STOPPED  (   ����           
   Data types                 Betriebs_Modus  *                  PID_Parameter  !   ����              Visualizations                myVisual  #   ����               Global Variables                 Global_Variables                     TwinCAT_Configuration  "                   Variable_Configuration  	   ����                                                             ��MR                         	   localhost            P      	   localhost            P      	   localhost            P            ��.�