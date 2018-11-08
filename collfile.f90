!       SUBROUTINE ORDER(NG,NTD,Y0,RAD,TIME)
! !
! !------------------------------------------------------------------------------
! ! Subroutine to list out the formation and destruction reactions for each 
! ! species in order of importance.
! !   IOTYP=0: Selected species (11) only
! !   IOTYP=1: All gas-phase species
! !   IOTYP=2: All species
! !------------------------------------------------------------------------------
!       IMPLICIT DOUBLE PRECISION(A-H,O-Z)
! !
!       PARAMETER(NSP=214,NCONS=2,NRM=2230,NSEL=16)
! !
!       CHARACTER*8 R1(NRM),R2(NRM),P1(NRM),P2(NRM),P3(NRM),P4(NRM), &
!      &            SPECI(NSP),SPSEL(NSEL),SPEC
!       REAL*8 G(NRM),R(NRM),X(NSP),Y0(NSP),XS(NCONS),TOT(NCONS),FLAG(10)
!       INTEGER INDXJ(NRM),LAB(10)
! !
!       COMMON/BLK1/SPECI
!       COMMON/BLK3/XS,TOT,R,D
!       COMMON/ORDER1/R1,R2,P1,P2,P3,P4
!       COMMON/ORDER2/INDXJ,IOTYP
! !
!       DATA  SPSEL/'N2H+','N2','NH3','ELECTR','CO','C','C+','HCO+', &
!      &      'CS','O','H2O','GNH3','GN2','GCO','GH2O','GO2'/
! !     DATA (SPSEL(I),I=1,33)/'H2+','H3+','CH','CH+','CH2','CH2+','CH3+',
! !    *     'CN','N2','N2H+','N+','NH','NH+','NH2','NH2+','NH3+','NH4+',
! !    *     'NH3','H2NC+','NO','HCO+','H3O+','OH','O2','CO','C2H+','H+',
! !    *     'HE+','C+','S+','C2S','HC2S+','ELECTR'/
!       DATA LOR/2/
! !
!       WRITE(2,4) RAD,TIME
!  4    FORMAT(1X,'Analysis at ',1PE8.2,' pc. and ',1PE8.2,' years',/)
! !
!       DO 112 K1=1,NRM
!          G(K1) = R(K1)
!  112  CONTINUE
!       DO 113 K1=1,NCONS
!          X(K1) = XS(K1)*D
!  113  CONTINUE
!       DO 114 K1=1,NTD
!          X(NCONS+K1) = Y0(K1)*D
!  114  CONTINUE
!       NTOT=NTD+NCONS
! !
!       DO 1 J=1,NRM
!         IF(INDXJ(J).EQ.9999) GO TO 2
!         DO 14 LI=1,NTOT
!            IF(R1(J).EQ.SPECI(LI)) G(J)=G(J)*X(LI)
!            IF(R2(J).EQ.SPECI(LI)) G(J)=G(J)*X(LI)
!  14     CONTINUE
!         IF(R1(J).EQ.'G') G(J)=G(J)*D
!         IF(R2(J).EQ.'G') G(J)=G(J)*D
!  1    CONTINUE
!  2    L=J-1
! !
! !------------------------------------------------------------------------------
! !
!       DO 90 I=1,L-1
!         DO 91 K=I+1,L
!           J1=INDXJ(I)
!           J2=INDXJ(K)
!           IF(G(J1).LT.G(J2)) THEN
!              INDXJ(I)=J2
!              INDXJ(K)=J1
!           END IF
!  91     CONTINUE
!  90   CONTINUE
! !
! !----Select species
! !
!       IF(IOTYP.EQ.0) LMAX=NSEL
!       IF(IOTYP.EQ.1) LMAX=NCONS+NG
!       IF(IOTYP.EQ.2) LMAX=NTOT
! !
!       DO 13 K=1,LMAX
!         IF(IOTYP.EQ.0) THEN
!            SPEC=SPSEL(K)
!            DO 18 ISEL=1,NTOT
!  18            IF(SPECI(ISEL).EQ.SPEC) GO TO 19
!  19        CONTINUE
!         ELSE 
!            SPEC=SPECI(K)
!            ISEL=K
!         END IF
!         FRATE = 0.0
!         DRATE = 0.0
!         K2=1
!         DO 16 K1=1,L
!           M=INDXJ(K1)
!           IF((R1(M).EQ.SPEC).OR.(R2(M).EQ.SPEC).OR.(P1(M).EQ.SPEC) .OR. &
!      &      (P2(M).EQ.SPEC).OR.(P3(M).EQ.SPEC).OR.(P4(M).EQ.SPEC)) THEN
!           FLAG(K2) = 1.0
!           IF((R1(M).EQ.SPEC).OR.(R2(M).EQ.SPEC)) FLAG(K2)= -1.0
!           IF(FLAG(K2).EQ.1.0) FRATE  = FRATE+G(M)
!           IF(FLAG(K2).EQ.-1.0) DRATE = DRATE-G(M)
!           LAB(K2)=M
!           K2=K2+1
!         END IF
!         IF(K2.GT.10) GO TO 17
!  16   CONTINUE
!  17   WRITE(LOR,9)
!       WRITE(LOR,11) SPEC,FRATE,-DRATE,X(ISEL)
!       WRITE(LOR,10)
!       DO 12 I=1,K2-1
!       J=LAB(I)
!       IR=0
!       IF((FLAG(I).EQ.1.0).AND.(FRATE.NE.0.0)) IR=NINT(G(J)*100.0/FRATE)
!       IF((FLAG(I).EQ.-1.0).AND.(DRATE.NE.0.0)) IR=NINT(G(J)*100.0/DRATE)
!       IF(IR.NE.0) WRITE(LOR,15) J,R1(J),R2(J),P1(J),P2(J),P3(J),P4(J),IR
!  12   CONTINUE
!  13   CONTINUE
! !
!       WRITE(LOR,20)
!  20   FORMAT(1X,78('-'),/)
!  15   FORMAT(1X,I4,1X,5(1A8,1X),A8,2X,I4,'%')
!  9    FORMAT(//,1X,5X,12('*'))
!  10   FORMAT(1X,5X,12('*'),/)
!  11   FORMAT(1X,5X,'*',1X,1A8,1X,'*',2X,'Frate=',1PE9.3,2X,'Drate=', &
!      &       1PE9.3,2X,'Abundance=',1PE9.3)
!       RETURN
!       END

! !******************************************************************************
! !******************************************************************************
! !******************************************************************************
! !
!       SUBROUTINE DIFFUN(N,T,Y,YDOT)

!         IMPLICIT DOUBLE PRECISION(A-H,O-Z)
!         REAL*8 PROD,LOSS
!   !
!         PARAMETER(NRM=2230,NSP=214,NCONS=2)
!   !
!         DIMENSION K(NRM),Y(N),YDOT(N)
!         COMMON /BLK3/X,TOTAL,K,D
!         COMMON /BLK9/FCO2,COCOV,FOX,PHOH,H2SHL,COSHL,FFRZ,FDES
!   !
!         CALL ADJUST(N,T,Y,YDOT)

!         INCLUDE 'rates/outputFiles/odes.f90'
!         ! LOSS = -RATE(10)*Y(11)*D-RATE(11)*Y(3)*D-RATE(12)*Y(35)&
!         !   &*D-RATE(13)*Y(102)*D-RATE(14)*Y(32)*D-RATE(112)*Y(119)*D-RATE(195)&
!         !   &*Y(54)*D-RATE(196)*Y(67)*D-RATE(197)*Y(4)*D-RATE(198)*Y(59)*D-RATE(199)&
!         !   &*Y(7)*D-RATE(200)*Y(27)*D-RATE(380)-RATE(423)-RATE(1117)*Y(183)&
!         !   &*D-RATE(1118)*Y(12)*D-RATE(1119)*Y(14)*D-RATE(1120)*Y(18)*D-RATE(1121)&
!         !   &*Y(23)*D-RATE(1122)*Y(29)*D-RATE(1123)*Y(117)*D-RATE(1124)*Y(120)&
!         !   &*D-RATE(1125)*Y(110)*D-RATE(1126)*Y(8)*D-RATE(1127)*Y(209)*D-RATE(1128)&
!         !   &*Y(81)*D-RATE(1129)*Y(196)*D-RATE(1757)*Y(41)*D-RATE(1758)*Y(51)&
!         !   &*D-RATE(1759)*Y(57)*D-RATE(1760)*Y(13)*D-RATE(1761)*Y(141)*D-RATE(1762)&
!         !   &*Y(17)*D-RATE(1763)*Y(22)*D-RATE(1764)*Y(11)*D-RATE(1765)*Y(150)&
!         !   &*D-RATE(1766)*Y(66)*D-RATE(1767)*Y(68)*D-RATE(1768)*Y(85)*D-RATE(1769)&
!         !   &*Y(35)*D-RATE(1770)*Y(116)*D-RATE(1771)*Y(58)*D-RATE(1772)*Y(76)&
!         !   &*D-RATE(1773)*Y(76)*D-RATE(1774)*Y(156)*D-RATE(1775)*Y(60)*D-RATE(1776)&
!         !   &*Y(92)*D-RATE(1777)*Y(92)*D-RATE(1778)*Y(92)*D-RATE(1779)*Y(109)&
!         !   &*D-RATE(1780)*Y(184)*D-RATE(1781)*Y(24)*D-RATE(1782)*Y(30)*D-RATE(1783)&
!         !   &*Y(19)*D-RATE(1784)*Y(168)*D-RATE(1785)*Y(87)*D-RATE(1786)*Y(87)&
!         !   &*D-RATE(1787)*Y(169)*D-RATE(1788)*Y(169)*D-RATE(1789)*Y(102)*D-RATE(1790&
!         !   &)*Y(111)*D-RATE(1791)*Y(111)*D-RATE(1792)*Y(111)*D-RATE(1793)&
!         !   &*Y(144)*D-RATE(1794)*Y(144)*D-RATE(1795)*Y(144)*D-RATE(1796)*Y(193)&
!         !   &*D-RATE(1797)*Y(32)*D-RATE(1798)*Y(207)*D-RATE(1799)*Y(175)*D-RATE(1800)&
!         !   &*Y(175)*D-RATE(1808)*Y(123)*D-RATE(2131)*Y(2)*D-RATE(2143)*Y(10)&
!         !   &*D-RATE(2144)*Y(9)*D-RATE(2145)*Y(26)*D-RATE(2146)*Y(32)*D-RATE(2147)&
!         !   &*Y(73)*D
!         ! PROD = +RATE(3)*Y(3)*Y(11)*D+2*RATE(4)*Y(3)*Y(3)*D+RATE(5)&
!         !   &*Y(3)*Y(35)*D+RATE(8)*Y(3)*Y(32)*D+2*RATE(9)*Y(3)*Y(215)*D+2*RATE(10)&
!         !   &*Y(1)*Y(11)*D+3*RATE(11)*Y(3)*Y(1)*D+2*RATE(12)*Y(1)*Y(35)*D+RATE(13)&
!         !   &*Y(1)*Y(102)*D+2*RATE(14)*Y(1)*Y(32)*D+RATE(113)*Y(2)*Y(118)*D+RATE(114)&
!         !   &*Y(2)*Y(41)*D+RATE(115)*Y(2)*Y(51)*D+RATE(116)*Y(2)*Y(47)*D+RATE(117)&
!         !   &*Y(2)*Y(126)*D+RATE(118)*Y(2)*Y(13)*D+RATE(119)*Y(2)*Y(17)*D+RATE(120)&
!         !   &*Y(2)*Y(22)*D+RATE(121)*Y(2)*Y(11)*D+RATE(122)*Y(151)*Y(2)*D+RATE(123)&
!         !   &*Y(2)*Y(85)*D+RATE(124)*Y(165)*Y(2)*D+RATE(125)*Y(2)*Y(35)*D+RATE(126)&
!         !   &*Y(2)*Y(116)*D+RATE(127)*Y(2)*Y(58)*D+RATE(128)*Y(2)*Y(76)*D+RATE(129)&
!         !   &*Y(2)*Y(123)*D+RATE(130)*Y(211)*Y(2)*D+RATE(131)*Y(2)*Y(109)*D+RATE(132)&
!         !   &*Y(2)*Y(43)*D+RATE(133)*Y(2)*Y(24)*D+RATE(134)*Y(2)*Y(30)*D+RATE(135)&
!         !   &*Y(19)*Y(2)*D+RATE(136)*Y(2)*Y(87)*D+RATE(137)*Y(2)*Y(169)*D+RATE(138)&
!         !   &*Y(2)*Y(102)*D+RATE(139)*Y(2)*Y(26)*D+RATE(140)*Y(2)*Y(193)*D+RATE(141)&
!         !   &*Y(2)*Y(32)*D+RATE(142)*Y(2)*Y(207)*D+RATE(143)*Y(2)*Y(104)*D+RATE(144)&
!         !   &*Y(2)*Y(210)*D+RATE(145)*Y(2)*Y(175)*D+RATE(146)*Y(2)*Y(72)*D+RATE(147)&
!         !   &*Y(2)*Y(185)*D+RATE(148)*Y(2)*Y(205)*D+RATE(149)*Y(133)*Y(2)*D+RATE(150)&
!         !   &*Y(2)*Y(89)*D+RATE(151)*Y(2)*Y(94)*D+RATE(152)*Y(106)*Y(2)*D+RATE(153)&
!         !   &*Y(2)*Y(80)*D+RATE(154)*Y(2)*Y(153)*D+RATE(155)*Y(195)*Y(2)*D+RATE(377)&
!         !   &*Y(3)+2*RATE(379)*Y(3)+RATE(386)*Y(51)+RATE(387)*Y(57)+RATE(391)&
!         !   &*Y(47)+RATE(398)*Y(12)+RATE(400)*Y(13)+RATE(402)*Y(17)+RATE(409)&
!         !   &*Y(11)+RATE(416)*Y(68)+RATE(419)*Y(35)+RATE(425)*Y(58)+RATE(426)&
!         !   &*Y(76)+RATE(430)*Y(123)+RATE(431)*Y(60)+RATE(433)*Y(92)+RATE(435)&
!         !   &*Y(109)+RATE(442)*Y(24)+RATE(443)*Y(30)+RATE(446)*Y(19)+RATE(455)&
!         !   &*Y(111)+RATE(460)*Y(32)+RATE(470)*Y(89)+RATE(471)*Y(94)+RATE(473)&
!         !   &*Y(80)+RATE(477)*Y(215)*Y(48)*D+2*RATE(479)*Y(215)*Y(52)*D+RATE(480)&
!         !   &*Y(215)*Y(52)*D+RATE(482)*Y(215)*Y(171)*D+RATE(484)*Y(215)*Y(171)&
!         !   &*D+RATE(485)*Y(215)*Y(171)*D+RATE(490)*Y(215)*Y(130)*D+RATE(492)&
!         !   &*Y(215)*Y(137)*D+RATE(495)*Y(215)*Y(12)*D+2*RATE(497)*Y(215)*Y(14)&
!         !   &*D+RATE(498)*Y(215)*Y(14)*D+RATE(499)*Y(215)*Y(18)*D+2*RATE(501)&
!         !   &*Y(215)*Y(18)*D+RATE(502)*Y(215)*Y(143)*D+RATE(504)*Y(215)*Y(108)&
!         !   &*D+RATE(506)*Y(215)*Y(108)*D+RATE(507)*Y(215)*Y(108)*D+RATE(508)&
!         !   &*Y(215)*Y(108)*D+2*RATE(509)*Y(215)*Y(23)*D+RATE(510)*Y(215)*Y(23)&
!         !   &*D+RATE(511)*Y(215)*Y(29)*D+2*RATE(513)*Y(215)*Y(29)*D+RATE(514)&
!         !   &*Y(215)*Y(29)*D+2*RATE(519)*Y(215)*Y(4)*D+2*RATE(522)*Y(215)*Y(86)&
!         !   &*D+RATE(523)*Y(215)*Y(86)*D+2*RATE(524)*Y(215)*Y(166)*D+RATE(525)&
!         !   &*Y(215)*Y(166)*D+2*RATE(526)*Y(215)*Y(125)*D+RATE(527)*Y(215)&
!         !   &*Y(125)*D+RATE(528)*Y(215)*Y(101)*D+2*RATE(531)*Y(215)*Y(36)*D+RATE(532)&
!         !   &*Y(215)*Y(36)*D+RATE(533)*Y(117)*Y(215)*D+2*RATE(534)*Y(117)*Y(215)&
!         !   &*D+RATE(535)*Y(215)*Y(214)*D+RATE(537)*Y(215)*Y(5)*D+3*RATE(538)&
!         !   &*Y(215)*Y(5)*D+RATE(541)*Y(215)*Y(91)*D+RATE(542)*Y(215)*Y(91)&
!         !   &*D+2*RATE(543)*Y(215)*Y(91)*D+RATE(544)*Y(215)*Y(172)*D+RATE(545)&
!         !   &*Y(215)*Y(172)*D+RATE(546)*Y(215)*Y(38)*D+RATE(547)*Y(215)*Y(38)&
!         !   &*D+2*RATE(549)*Y(215)*Y(38)*D+RATE(550)*Y(215)*Y(120)*D+2*RATE(552)&
!         !   &*Y(215)*Y(120)*D+RATE(553)*Y(215)*Y(120)*D+RATE(555)*Y(215)*Y(197)&
!         !   &*D+RATE(556)*Y(215)*Y(59)*D+2*RATE(557)*Y(215)*Y(69)*D+RATE(558)&
!         !   &*Y(215)*Y(69)*D+RATE(559)*Y(215)*Y(69)*D+RATE(560)*Y(215)*Y(77)&
!         !   &*D+RATE(561)*Y(215)*Y(155)*D+RATE(562)*Y(215)*Y(155)*D+RATE(565)&
!         !   &*Y(215)*Y(157)*D+RATE(566)*Y(124)*Y(215)*D+RATE(567)*Y(93)*Y(215)&
!         !   &*D+RATE(568)*Y(215)*Y(173)*D+RATE(569)*Y(215)*Y(78)*D+RATE(571)&
!         !   &*Y(215)*Y(198)*D+RATE(572)*Y(215)*Y(110)*D+RATE(574)*Y(215)*Y(212)&
!         !   &*D+RATE(575)*Y(215)*Y(178)*D+RATE(576)*Y(215)*Y(213)*D+RATE(577)&
!         !   &*Y(215)*Y(213)*D+RATE(580)*Y(215)*Y(199)*D+RATE(581)*Y(215)*Y(8)&
!         !   &*D+RATE(583)*Y(215)*Y(79)*D+RATE(585)*Y(215)*Y(20)*D+2*RATE(586)&
!         !   &*Y(215)*Y(25)*D+RATE(587)*Y(215)*Y(25)*D+RATE(588)*Y(31)*Y(215)&
!         !   &*D+2*RATE(589)*Y(31)*Y(215)*D+2*RATE(591)*Y(215)*Y(37)*D+RATE(592)&
!         !   &*Y(215)*Y(37)*D+RATE(596)*Y(112)*Y(215)*D+RATE(600)*Y(215)*Y(33)&
!         !   &*D+RATE(610)*Y(215)*Y(81)*D+2*RATE(612)*Y(90)*Y(215)*D+RATE(613)&
!         !   &*Y(90)*Y(215)*D+RATE(614)*Y(215)*Y(95)*D+RATE(617)*Y(215)*Y(107)&
!         !   &*D+RATE(619)*Y(215)*Y(113)*D+RATE(622)*Y(215)*Y(158)*D+RATE(625)&
!         !   &*Y(10)*Y(47)*D+RATE(626)*Y(10)*Y(13)*D+RATE(628)*Y(10)*Y(17)*D+RATE(633)&
!         !   &*Y(10)*Y(11)*D+RATE(638)*Y(10)*Y(35)*D+RATE(639)*Y(10)*Y(35)*D+RATE(640)&
!         !   &*Y(10)*Y(116)*D+RATE(643)*Y(10)*Y(182)*D+RATE(645)*Y(10)*Y(60)&
!         !   &*D+RATE(646)*Y(10)*Y(109)*D+RATE(647)*Y(10)*Y(24)*D+RATE(649)&
!         !   &*Y(10)*Y(19)*D+RATE(655)*Y(10)*Y(32)*D+RATE(662)*Y(10)*Y(80)*D+RATE(702)&
!         !   &*Y(73)*Y(47)*D+RATE(703)*Y(9)*Y(48)*D+RATE(704)*Y(9)*Y(12)*D+RATE(705)&
!         !   &*Y(9)*Y(14)*D+RATE(709)*Y(117)*Y(9)*D+RATE(715)*Y(110)*Y(9)*D+RATE(721)&
!         !   &*Y(9)*Y(81)*D+RATE(723)*Y(41)*Y(12)*D+RATE(729)*Y(12)*Y(22)*D+RATE(731)&
!         !   &*Y(12)*Y(53)*D+RATE(736)*Y(35)*Y(12)*D+RATE(742)*Y(12)*Y(58)*D+RATE(746)&
!         !   &*Y(12)*Y(15)*D+RATE(753)*Y(12)*Y(26)*D+RATE(757)*Y(104)*Y(12)&
!         !   &*D+RATE(761)*Y(35)*Y(14)*D+RATE(762)*Y(14)*Y(116)*D+RATE(764)&
!         !   &*Y(14)*Y(116)*D+RATE(768)*Y(14)*Y(26)*D+RATE(771)*Y(104)*Y(14)&
!         !   &*D+RATE(790)*Y(105)*Y(13)*D+RATE(801)*Y(26)*Y(18)*D+RATE(808)&
!         !   &*Y(105)*Y(17)*D+RATE(824)*Y(22)*Y(52)*D+RATE(834)*Y(71)*Y(22)&
!         !   &*D+RATE(839)*Y(105)*Y(22)*D+RATE(840)*Y(105)*Y(22)*D+RATE(852)&
!         !   &*Y(43)*Y(29)*D+RATE(855)*Y(42)*Y(11)*D+RATE(870)*Y(16)*Y(11)*D+RATE(875)&
!         !   &*Y(11)*Y(27)*D+RATE(879)*Y(105)*Y(11)*D+RATE(880)*Y(73)*Y(11)&
!         !   &*D+RATE(884)*Y(58)*Y(54)*D+RATE(901)*Y(2)*Y(65)*D+RATE(910)*Y(2)&
!         !   &*Y(85)*D+RATE(913)*Y(2)*Y(116)*D+RATE(927)*Y(41)*Y(4)*D+RATE(929)&
!         !   &*Y(47)*Y(4)*D+RATE(930)*Y(9)*Y(4)*D+RATE(931)*Y(13)*Y(4)*D+RATE(932)&
!         !   &*Y(22)*Y(4)*D+RATE(933)*Y(22)*Y(4)*D+RATE(934)*Y(11)*Y(4)*D+RATE(935)&
!         !   &*Y(53)*Y(4)*D+RATE(936)*Y(150)*Y(4)*D+RATE(937)*Y(66)*Y(4)*D+RATE(938)&
!         !   &*Y(3)*Y(4)*D+RATE(939)*Y(85)*Y(4)*D+RATE(940)*Y(35)*Y(4)*D+RATE(941)&
!         !   &*Y(116)*Y(4)*D+RATE(944)*Y(4)*Y(6)*D+RATE(945)*Y(70)*Y(4)*D+RATE(946)&
!         !   &*Y(4)*Y(15)*D+RATE(947)*Y(19)*Y(4)*D+RATE(948)*Y(4)*Y(87)*D+RATE(949)&
!         !   &*Y(102)*Y(4)*D+RATE(950)*Y(26)*Y(4)*D+RATE(951)*Y(4)*Y(32)*D+RATE(952)&
!         !   &*Y(3)*Y(10)*D+RATE(953)*Y(3)*Y(42)*D+RATE(954)*Y(3)*Y(48)*D+RATE(955)&
!         !   &*Y(3)*Y(12)*D+RATE(956)*Y(3)*Y(14)*D+RATE(957)*Y(3)*Y(23)*D+RATE(958)&
!         !   &*Y(3)*Y(54)*D+RATE(959)*Y(3)*Y(67)*D+RATE(960)*Y(3)*Y(67)*D+RATE(961)&
!         !   &*Y(3)*Y(152)*D+RATE(962)*Y(3)*Y(119)*D+RATE(963)*Y(3)*Y(36)*D+RATE(964)&
!         !   &*Y(3)*Y(117)*D+RATE(965)*Y(3)*Y(59)*D+RATE(966)*Y(3)*Y(124)*D+RATE(967)&
!         !   &*Y(3)*Y(110)*D+RATE(968)*Y(3)*Y(7)*D+RATE(970)*Y(3)*Y(16)*D+RATE(971)&
!         !   &*Y(3)*Y(71)*D+RATE(973)*Y(3)*Y(20)*D+RATE(974)*Y(3)*Y(25)*D+RATE(975)&
!         !   &*Y(3)*Y(31)*D+RATE(976)*Y(3)*Y(27)*D+RATE(978)*Y(3)*Y(33)*D+RATE(979)&
!         !   &*Y(3)*Y(105)*D+RATE(980)*Y(3)*Y(209)*D+RATE(981)*Y(3)*Y(107)*D+RATE(982)&
!         !   &*Y(3)*Y(154)*D+RATE(990)*Y(85)*Y(103)*D+RATE(1006)*Y(104)*Y(36)&
!         !   &*D+RATE(1031)*Y(73)*Y(35)*D+RATE(1072)*Y(5)*Y(43)*D+RATE(1080)&
!         !   &*Y(5)*Y(39)*D+RATE(1082)*Y(5)*Y(26)*D+RATE(1111)*Y(39)*Y(38)*D+RATE(1139&
!         !   &)*Y(58)*Y(122)*D+RATE(1199)*Y(51)*Y(7)*D+RATE(1202)*Y(57)*Y(7)&
!         !   &*D+RATE(1203)*Y(65)*Y(7)*D+RATE(1206)*Y(47)*Y(7)*D+RATE(1213)&
!         !   &*Y(13)*Y(7)*D+RATE(1222)*Y(22)*Y(7)*D+RATE(1224)*Y(22)*Y(7)*D+RATE(1226)&
!         !   &*Y(11)*Y(7)*D+RATE(1237)*Y(85)*Y(7)*D+RATE(1242)*Y(35)*Y(7)*D+RATE(1244)&
!         !   &*Y(116)*Y(7)*D+RATE(1246)*Y(167)*Y(7)*D+RATE(1249)*Y(58)*Y(7)&
!         !   &*D+RATE(1251)*Y(58)*Y(7)*D+RATE(1253)*Y(76)*Y(7)*D+RATE(1257)&
!         !   &*Y(156)*Y(7)*D+RATE(1259)*Y(123)*Y(7)*D+RATE(1260)*Y(60)*Y(7)&
!         !   &*D+RATE(1261)*Y(60)*Y(7)*D+RATE(1263)*Y(92)*Y(7)*D+RATE(1266)&
!         !   &*Y(211)*Y(7)*D+RATE(1267)*Y(109)*Y(7)*D+RATE(1271)*Y(24)*Y(7)&
!         !   &*D+RATE(1273)*Y(30)*Y(7)*D+RATE(1274)*Y(19)*Y(7)*D+RATE(1286)&
!         !   &*Y(32)*Y(7)*D+RATE(1297)*Y(89)*Y(7)*D+RATE(1299)*Y(94)*Y(7)*D+RATE(1301)&
!         !   &*Y(106)*Y(7)*D+RATE(1302)*Y(80)*Y(7)*D+RATE(1307)*Y(100)*Y(16)&
!         !   &*D+RATE(1309)*Y(100)*Y(16)*D+RATE(1310)*Y(100)*Y(16)*D+RATE(1311)&
!         !   &*Y(16)*Y(22)*D+RATE(1312)*Y(16)*Y(22)*D+2*RATE(1313)*Y(16)*Y(22)&
!         !   &*D+RATE(1320)*Y(16)*Y(116)*D+RATE(1326)*Y(19)*Y(16)*D+RATE(1332)&
!         !   &*Y(71)*Y(85)*D+RATE(1333)*Y(71)*Y(116)*D+RATE(1344)*Y(48)*Y(15)&
!         !   &*D+RATE(1347)*Y(52)*Y(15)*D+RATE(1349)*Y(14)*Y(15)*D+RATE(1351)&
!         !   &*Y(36)*Y(15)*D+RATE(1354)*Y(110)*Y(15)*D+RATE(1355)*Y(20)*Y(15)&
!         !   &*D+RATE(1356)*Y(25)*Y(15)*D+RATE(1358)*Y(33)*Y(15)*D+RATE(1364)&
!         !   &*Y(41)*Y(20)*D+RATE(1391)*Y(104)*Y(20)*D+RATE(1410)*Y(104)*Y(25)&
!         !   &*D+RATE(1412)*Y(42)*Y(24)*D+RATE(1463)*Y(19)*Y(42)*D+RATE(1475)&
!         !   &*Y(19)*Y(27)*D+RATE(1479)*Y(19)*Y(105)*D+RATE(1481)*Y(39)*Y(91)&
!         !   &*D+RATE(1482)*Y(39)*Y(69)*D+RATE(1483)*Y(39)*Y(69)*D+RATE(1501)&
!         !   &*Y(27)*Y(32)*D+RATE(1503)*Y(103)*Y(51)*D+RATE(1504)*Y(100)*Y(103)&
!         !   &*D+RATE(1522)*Y(157)*Y(26)*D+RATE(1525)*Y(110)*Y(26)*D+RATE(1528)&
!         !   &*Y(25)*Y(26)*D+RATE(1532)*Y(26)*Y(33)*D+RATE(1534)*Y(26)*Y(81)&
!         !   &*D+RATE(1535)*Y(90)*Y(26)*D+RATE(1555)*Y(104)*Y(33)*D+RATE(1564)&
!         !   &*Y(77)*Y(32)*D+RATE(1569)*Y(105)*Y(32)*D+RATE(1570)*Y(73)*Y(32)&
!         !   &*D+RATE(1571)*Y(105)*Y(116)*D+RATE(1574)*Y(120)*Y(104)*D+RATE(1589)&
!         !   &*Y(90)*Y(104)*D+RATE(1590)*Y(105)*Y(80)*D+RATE(1591)*Y(41)*Y(51)&
!         !   &*D+RATE(1592)*Y(41)*Y(58)*D+RATE(1595)*Y(51)*Y(87)*D+RATE(1599)&
!         !   &*Y(47)*Y(58)*D+RATE(1600)*Y(47)*Y(60)*D+RATE(1603)*Y(9)*Y(57)&
!         !   &*D+RATE(1604)*Y(9)*Y(75)*D+RATE(1606)*Y(9)*Y(128)*D+RATE(1607)&
!         !   &*Y(9)*Y(13)*D+RATE(1609)*Y(9)*Y(17)*D+RATE(1610)*Y(9)*Y(11)*D+RATE(1616)&
!         !   &*Y(109)*Y(9)*D+RATE(1620)*Y(9)*Y(24)*D+RATE(1621)*Y(9)*Y(24)*D+RATE(1623&
!         !   &)*Y(19)*Y(9)*D+RATE(1632)*Y(9)*Y(32)*D+RATE(1638)*Y(9)*Y(80)*D+2&
!         !   &*RATE(1640)*Y(13)*Y(13)*D+RATE(1641)*Y(13)*Y(13)*D+RATE(1652)&
!         !   &*Y(13)*Y(87)*D+2*RATE(1654)*Y(13)*Y(102)*D+2*RATE(1659)*Y(13)&
!         !   &*Y(26)*D+RATE(1660)*Y(13)*Y(26)*D+RATE(1662)*Y(13)*Y(32)*D+RATE(1666)&
!         !   &*Y(13)*Y(104)*D+RATE(1669)*Y(17)*Y(17)*D+RATE(1685)*Y(17)*Y(26)&
!         !   &*D+RATE(1686)*Y(17)*Y(26)*D+RATE(1690)*Y(104)*Y(17)*D+RATE(1695)&
!         !   &*Y(11)*Y(51)*D+RATE(1696)*Y(11)*Y(65)*D+RATE(1697)*Y(11)*Y(22)&
!         !   &*D+RATE(1703)*Y(11)*Y(15)*D+RATE(1707)*Y(11)*Y(87)*D+RATE(1708)&
!         !   &*Y(11)*Y(102)*D+RATE(1709)*Y(11)*Y(102)*D+RATE(1714)*Y(11)*Y(26)&
!         !   &*D+RATE(1716)*Y(193)*Y(11)*D+RATE(1717)*Y(11)*Y(32)*D+RATE(1718)&
!         !   &*Y(104)*Y(11)*D+RATE(1721)*Y(11)*Y(175)*D+RATE(1722)*Y(53)*Y(51)&
!         !   &*D+RATE(1726)*Y(58)*Y(53)*D+RATE(1728)*Y(60)*Y(53)*D+RATE(1741)&
!         !   &*Y(3)*Y(118)*D+RATE(1742)*Y(3)*Y(47)*D+RATE(1743)*Y(3)*Y(9)*D+RATE(1744)&
!         !   &*Y(3)*Y(13)*D+RATE(1745)*Y(3)*Y(17)*D+RATE(1746)*Y(3)*Y(11)*D+RATE(1747)&
!         !   &*Y(3)*Y(53)*D+RATE(1748)*Y(3)*Y(109)*D+RATE(1749)*Y(3)*Y(15)*D+RATE(1750&
!         !   &)*Y(3)*Y(24)*D+RATE(1751)*Y(3)*Y(19)*D+RATE(1752)*Y(3)*Y(102)&
!         !   &*D+RATE(1754)*Y(3)*Y(26)*D+RATE(1755)*Y(3)*Y(32)*D+RATE(1756)&
!         !   &*Y(3)*Y(104)*D+RATE(1775)*Y(1)*Y(60)*D+RATE(1816)*Y(47)*Y(15)&
!         !   &*D+RATE(1818)*Y(128)*Y(15)*D+RATE(1820)*Y(177)*Y(15)*D+RATE(1822)&
!         !   &*Y(13)*Y(15)*D+RATE(1823)*Y(13)*Y(15)*D+RATE(1825)*Y(17)*Y(15)&
!         !   &*D+2*RATE(1827)*Y(17)*Y(15)*D+RATE(1834)*Y(76)*Y(15)*D+RATE(1837)&
!         !   &*Y(109)*Y(15)*D+RATE(1840)*Y(19)*Y(15)*D+RATE(1848)*Y(32)*Y(15)&
!         !   &*D+RATE(1856)*Y(87)*Y(24)*D+2*RATE(1865)*Y(19)*Y(19)*D+RATE(1868)&
!         !   &*Y(19)*Y(87)*D+RATE(1872)*Y(19)*Y(26)*D+RATE(1875)*Y(19)*Y(32)&
!         !   &*D+RATE(1878)*Y(19)*Y(104)*D+RATE(1890)*Y(57)*Y(26)*D+RATE(1912)&
!         !   &*Y(58)*Y(26)*D+RATE(1913)*Y(76)*Y(26)*D+RATE(1916)*Y(156)*Y(26)&
!         !   &*D+RATE(1917)*Y(92)*Y(26)*D+RATE(1921)*Y(109)*Y(26)*D+RATE(1923)&
!         !   &*Y(26)*Y(24)*D+RATE(1935)*Y(26)*Y(32)*D+2*RATE(1944)*Y(89)*Y(26)&
!         !   &*D+RATE(1945)*Y(26)*Y(94)*D+RATE(1947)*Y(80)*Y(26)*D+RATE(1949)&
!         !   &*Y(51)*Y(32)*D+RATE(1954)*Y(53)*Y(32)*D+RATE(1955)*Y(66)*Y(32)&
!         !   &*D+RATE(1957)*Y(151)*Y(32)*D+RATE(1966)*Y(87)*Y(32)*D+RATE(1969)&
!         !   &*Y(104)*Y(32)*D+RATE(1970)*Y(175)*Y(32)*D+RATE(1971)*Y(72)*Y(32)&
!         !   &*D+RATE(1972)*Y(104)*Y(76)*D+RATE(1974)*Y(109)*Y(104)*D+RATE(1984)&
!         !   &*Y(48)+RATE(1986)*Y(51)+RATE(1987)*Y(57)+RATE(1991)*Y(47)+RATE(2000)&
!         !   &*Y(14)+RATE(2003)*Y(13)+RATE(2006)*Y(18)+RATE(2007)*Y(17)+RATE(2012)&
!         !   &*Y(100)+RATE(2015)*Y(23)+RATE(2017)*Y(22)+RATE(2019)*Y(22)+RATE(2020)&
!         !   &*Y(11)+RATE(2030)*Y(4)+RATE(2031)*Y(68)+2*RATE(2033)*Y(85)+RATE(2035)&
!         !   &*Y(85)+RATE(2037)*Y(36)+RATE(2039)*Y(35)+RATE(2041)*Y(116)+2*RATE(2044)&
!         !   &*Y(167)+RATE(2045)*Y(5)+RATE(2048)*Y(58)+RATE(2049)*Y(77)+RATE(2050)&
!         !   &*Y(76)+RATE(2054)*Y(123)+RATE(2056)*Y(60)+RATE(2058)*Y(92)+RATE(2059)&
!         !   &*Y(110)+RATE(2063)*Y(109)+RATE(2070)*Y(24)+RATE(2071)*Y(30)+RATE(2074)&
!         !   &*Y(19)+RATE(2084)*Y(111)+RATE(2089)*Y(33)+RATE(2090)*Y(32)+RATE(2103)&
!         !   &*Y(81)+RATE(2105)*Y(89)+RATE(2106)*Y(94)+RATE(2110)*Y(106)+RATE(2111)&
!         !   &*Y(106)+RATE(2112)*Y(80)+RATE(2156)*Y(215)*Y(2)*D+RATE(2170)*Y(137)&
!         !   &*D+RATE(2173)*Y(171)*D+RATE(2174)*Y(108)*D+RATE(2189)*Y(197)*D+RATE(2213&
!         !   &)*Y(199)*D+RATE(2215)*Y(113)*D+RATE(2220)*Y(125)*D+RATE(2264)&
!         !   &*Y(38)*D+RATE(2265)*Y(155)*D+RATE(2266)*Y(29)*D+RATE(2291)*Y(101)&
!         !   &*D+RATE(2296)*Y(120)*D+RATE(2298)*Y(172)*D+RATE(2299)*Y(178)*D+RATE(2301&
!         !   &)*Y(198)*D+RATE(2304)*Y(37)*D+RATE(2305)*Y(69)*D+RATE(2306)*Y(79)&
!         !   &*D+RATE(2307)*Y(173)*D+RATE(2311)*Y(213)*D
!         ! YDOT(1) = PROD+Y(1)*LOSS

!         RETURN
!       END
! !