function batcher
%
% 3mm/s oro-anal propagation; 
% 100um ICCmy compartment;
% 33ms ICCmy compartment transit time;
% 631um is small intestine LM length constant
%

clear global

global DURN STARTS LAMBDA GAPS PULON PULOFF IPULSE VM SM VH SH
global LAMINC TAUH1 TAUHLIM TAUH2 NUMUNITS UNITS DETAILS NEWNOISE
global KF KB KM KH KL STH FIGNUM GPRIMAX TAUACT1 VACT SACT VINACT SINACT
global MBAS HBAS VTOP VBOT AU BU GPLATM AMP RANDAMPS HILL RICC AX1
global COMBAS COHBAS TBOT LGAPS LSTARTS
global TAUIN1 STIN PROSPECTS MAXINTSTEP RATEFAC UNITSCAL
global MAXPROS PREV IDENTAMPS AMPINC
global NUMCOM PASSIVE IPULSEL PULONL PULOFFL PULONL2 PULOFFL2
global ROUND TAUINLIM TAUIN2 TAUM1 NUDGE KFP KBP LBOT LTOP
global PEELED TRANS MAXXO MAXY STANDICS

pack;
tic;

NUMCOM = 20;
PASSIVE = 1; %1; %0;
ROUND = 1;
PEELED = 0;
TRANS = 0;   %0.1       %0,1=not transverse; 0<TRANS<1 = transverse gLMLM
DETAILS = 0;
DURN = 50; %50; %120;              %31;      % 1200 works ok on Dell 3GHz;
FIGNUM = 50;
MAXINTSTEP = 3e-3;      %6e-3; %3e-3; %6e-4 differs from 9e-4;  % sec
PROSPECTS = 20;         %18; %30; %18;
IDENTAMPS = 0;          % can be set when PROSPECTS or DURN is varied
STANDICS = 0; %1;       % -1 reads/writes FCs; 0 reads FCs; 1 uses equal values; 2 uses equal values and writes FCs;
VTOP = 0;
VBOT = -80;
TBOT = 0;
LBOT = 0.5;
LTOP = 0.9;
%
%   Unit parameters
%
UNITS = 1;             %1; %0;
RANDAMPS = 1;
NEWNOISE = 0;
RICC = 12.3e-3;         %27.2e-3 for small bundle;      % gigohms  see Cousins et al 2003, Table 1
%       12.3e-3 gives ~100uV minimum unit in CMB, ~300uV in ICCmy
UNITSCAL = 1;           %1;       % multiplies unit size, RICC independent
RATEFAC = 2;            %2;       % 1 for CMB; multiplies max frequency
if (UNITS == 0)
    RATEFAC = 6 * RATEFAC;
end
        
unitbase = 0.2;          %0.2;          % 0.2,25 gives 0.5mV unit
unscal = UNITSCAL * unitbase;
GPLATM = unscal / RICC;  % multiplies unit size (see  above)
AU = -1.0 / 0.434;       % average of 17 Edw, Hir & Suz
BU = -1.0 / 0.077;       % average of 17 Edw, Hir & Suz
%
%   primary component voltage dependencies & kinetics
%
vshift = 0; %-6; %-6;
NUDGE = -10; %0;            % -10 for drive; 0 for no drive
GPRIMAX = 1.2*2.035e7;      %1.2*2.035e7    %0.94* slows rate, not vel much

KFP = 0.012;            %0.012;    %0.0114 slows velocity & rate;
KBP = 6;                %6;        %6.5 speeds velocity and slows rate

TAUACT1 = 0.06;         %0.06;     % 0.06 gives pointier ears than 0.1; ok for UNITS=1
TAUIN1 = 0.25;          %0.25;
TAUINLIM = -48 + vshift;  %-48 + vshift;
STIN = 2;               %2;
TAUIN2 = 90;            %90;            % 3 for Hiroe's spikes during diastole
VACT = -48 + vshift;    %-48 + vshift;
SACT = 0.78;            %0.78;           % 0.9 slows rate, 0.7 speeds rate
VINACT = -53 + vshift;  %-53 + vshift;   % -50 speeds rate, -56 slows rate 
SINACT = 1;             %1;              % 2 speeds rate
%
%   messenger formation reactions
%
lammax = RATEFAC * 28 / unscal;
LAMINC = lammax;
AMPINC = 9;

KF = 5.5;               % 5.5;
KF = KF * 0.0157;
KB = 3.3;               % 3.3;
KM = 0.4478;            % 0.44776;
KH = 0.3;               % 0.3;
KL = 7.5e-4;            % 7.5e-4;
HILL = 6.8;             % 6.8;
%
%   KF voltage dependencies & kinetics
%
MBAS = 0.0;
VM = -59.5;                 %-59.5;
SM = 0.5;                   %0.5;   % 0.3 for CMB
HBAS = 0.0;
VH = -67.5;                 %-67.5;
SH = 0.58;                  %0.58;  % 0.32 for CMB;

TAUM1 = 0.15;               %0.15 ok for UNITS=1;  %used to be 0.25; see earlier paper
TAUH1 = 1.7;                %1.7; %2.1; %2.5; %3; 
TAUHLIM = -45.5;            %-45.5
STH = 0.5;                  %0.5
TAUH2 = 6;                  %6

%   current injection into ICCmy

IPULSE = 0;                 %-1000; %7000; %0;
PULON = 3;                  %33; %5; %9999;
PULOFF = 999;               %33.1; %5.1; %10000;

%   current injections into LM

IPULSEL = 0; %80000;        %1000 ; %2000; %0;
PULONL = 10;                %5; %9999;
PULOFFL = 10.05;            %5.1; %10000;
PULONL2 = 57;               %9999;
PULOFFL2 = 58;              %10000;
 
COMBAS = 1 - MBAS;
COHBAS = 1 - HBAS;
NUMUNITS = 0;
MAXXO = 0;
MAXY = 0;
for ii = 1:NUMCOM
    MAXPROS(ii) = 0;
    PREV(ii) = 0;
end

if (PASSIVE > 0)
    IPULSE = 0;
    if (PEELED > 0)
        IPULSE = 10000;
    end
    PULON = 5;
    PULOFF = 10;
    IPULSEL = 10000;        %0;
    PULONL = 5;
    PULOFFL = 10;
    DURN = 10;
    GPRIMAX = 0;
    KF = 0;
    STANDICS = 1;           % 1 uses equal values;
end

if (NEWNOISE == 0)
   rand('state', 0);
end

d1i;

if (NEWNOISE ~= 0)
    NEWNOISE
end

beep on
beep
beep off

% the batching loop went here

% --------------------------------------------------------------------------

function d1i
%
%   one ICC, lumped Rm, no ions, Hodgkin-Huxley kF
%
%                   kF(Em,t)         kM     ^ kH   kL(sigmoid threshold)
%      [Precursors] <> [Intermediate] > [Messenger] > Unit Rate > gPlat(Em,t)
%                   kB
%
%      gPlat(Em,t) & gPrim(Em,t) > Em 
%
%   Units:  mV, gOhm, nS, pA, nF, sec, mM, uL
%
global DURN STARTS LAMBDA GAPS IPULSE LAMINC NUMUNITS FIGNUM UNITS DETAILS
global VTOP VBOT TBOT LGAPS LSTARTS RANDAMPS PROSPECTS MAXINTSTEP MAXPROS
global NUMCOM PASSIVE HILL KL RATEFAC IPULSEL LBOT LTOP MAXXO MAXY STANDICS

vi0 = -65;         % -65;
m0 = 0.04;         % 0.04;
h0 = 0.001;        % 0.001;
interm0 = 40e-6;   % 40e-6;
mess0 = 200e-6;    % 200e-6;
act0 = 0.01;       % 0.01;        % 0.0033 for CMB;
inact0 = 0.001;    % 0.001;       % 0.2489 for CMB;
messP0 = 1e-7;     % 1e-7;
vLong0 = -65;      % -65;

LAMBDA = 0.3;      % 0.3;
evlen = 0;

if (DETAILS == 1)
    figure(FIGNUM);
end

    if (STANDICS > 0)
        y0 = [vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0...
            ];
    else
%        fid = fopen('C:\Matfiles\ICs180.dat');
        fid = fopen('ICs180.dat');
        y0 = fscanf(fid,'%g',[180 inf]);    % It has 180 rows now.
        y0 = y0';
        fclose(fid);
    end

if (UNITS > 0)
    if (DURN <= 20)
        evlen = round(LAMINC * DURN + PROSPECTS);
    else
        evlen = round(0.6 * LAMINC * DURN + PROSPECTS);   % 0.6
    end
    for ii = 1:NUMCOM
        GAPS(ii,:) = -log(rand(1,evlen));
        LGAPS(ii,:) = length(GAPS(ii,:));
        STARTS(ii,:) = cumsum(GAPS(ii,:) / LAMBDA);
        LSTARTS(ii,:) = length(STARTS(ii,:));
    end

    if (RANDAMPS > 0)
        getamps(evlen);
    end
end

    options = [];
    options = odeset('MaxStep', MAXINTSTEP);
    [t,y] = ode15s(@f,[0 DURN],y0,options);

% inacts;
if (DETAILS == 1)
% this plots m & h
    figure(210);
    plot(t,y(:,2),t,y(:,3));
%    subplot(6,1,4); plot(t,y(:,2),t,y(:,3),t,y(:,2).*y(:,3));
    axis([TBOT DURN 0 1]);
    intp = 10 .* y(:,4);
    messp = 10 .* y(:,5);
% this plots [Intermediate] & [Messenger]
    figure(220);
    plot(t,y(:,4),t,y(:,5));
%    subplot(6,1,5); plot(t,y(:,4),t,y(:,5),t,intp,t,messp);
    axis([TBOT DURN 0 0.005]);
% this plots LAMBDA    
    ttemp = y(:,5) .^ HILL;
    ttemp2 = KL ^ HILL;
    urate = (LAMINC * ttemp) ./ (ttemp2 + ttemp);
    figure(240);
    plot(t,urate);
    axis([TBOT DURN 0 150*RATEFAC]);
% end
% this plots activation & inactivation for Primary component
    figure(270);
    plot(t,y(:,6),t,y(:,7),t,y(:,6) .* y(:,7));
   axis([TBOT DURN -0.05 1.05]);
end

if (PASSIVE > 0)
%    figure(FIGNUM);
%    plot(t,y(:,1),t,y(:,10),t,y(:,19),t,y(:,28),t,y(:,37),t,y(:,46),t,y(:,172));
%    axis([TBOT DURN VBOT VTOP]);
%    figure(FIGNUM+1);
%    plot(t,y(:,9),t,y(:,18),t,y(:,27),t,y(:,36),t,y(:,45),t,y(:,54),t,y(:,180));
%    axis([TBOT DURN VBOT VTOP]);
else
   figure(FIGNUM);
   plot(t,y(:,1),t,y(:,10),t,y(:,19),t,y(:,28),t,y(:,37),t,y(:,46),t,y(:,172));
   axis([TBOT DURN VBOT VTOP]);
   figure(FIGNUM+1);
%    plot(t,y(:,1),t,y(:,19),t,y(:,37),t,y(:,55),t,y(:,73),t,y(:,91),t,y(:,109),t,y(:,127),t,y(:,145),t,y(:,163),t,y(:,172));
%    axis([DURN-20 DURN -45 -20]);
%    figure(FIGNUM+2);
   plot(t,y(:,19),t,y(:,109));
   axis([TBOT DURN VBOT VTOP]);
%    axis([DURN-20 DURN -45 -20]);
%    axis([0 DURN VBOT VTOP]);
%    plot(t,y(:,9),t,y(:,18),t,y(:,27),t,y(:,36),t,y(:,45),t,y(:,54),t,y(:,180));
%    axis([30 DURN -65 -34]);
% Plots y(8) [messenger] for GPrim
%    figure(FIGNUM+2);
%    plot(t,y(:,8),t,y(:,44),t,y(:,80),t,y(:,116),t,y(:,152),t,y(:,179));
%    axis([0 DURN 0 2e-4]);
end
    
%    plot(t,y(:,9),t,y(:,18),t,y(:,27),t,y(:,36),t,y(:,45),t,y(:,54));

if (PASSIVE > 0)
    len = length(y(:,1));
    prepwidth = 0.1;
    preplength = 0.6;
    if (IPULSE~=0 | IPULSEL~=0)
        for i = 1:NUMCOM-1
            over = 1 + 9*i;
            under = over - 9;
            lcMY(i) = prepwidth / (-log((y(len-1,over) - y(2,over)) / (y(len-1,under) - y(2,under))));
            lcMYone(i) = prepwidth * i / (-log((y(len-1,over) - y(2,over)) / (y(len-1,1) - y(2,1))));
            over = over + 8;
            under = over - 9;
            lcLM(i) = prepwidth / (-log((y(len-1,over) - y(2,over)) / (y(len-1,under) - y(2,under))));
            lcLMone(i) = prepwidth * i / (-log((y(len-1,over) - y(2,over)) / (y(len-1,9) - y(2,9))));
            pos(i) = i;
        end
        figure(FIGNUM);
        plot(pos,lcLM,pos,lcLMone,pos,lcMY,pos,lcMYone);
        grid on
        axis([0 NUMCOM LBOT LTOP]);
    end
end

if (NUMUNITS > evlen)
    NUMUNITS
    evlen
end

for ii = 1:NUMCOM
    if (MAXPROS(ii) > PROSPECTS)
        ii
        MAXPROS(ii)
        PROSPECTS
    end    
end

if ((STANDICS < 0) | (STANDICS == 2))
%    fid = fopen('C:\Matfiles\ICs180.dat','w');
    fid = fopen('ICs180.dat','w');
    fprintf(fid,'%30.25f\n',y(end,:));
    fclose(fid);
end

% fid = fopen('C:\Matfiles\FCs180.dat','w');
fid = fopen('FCs180.dat','w');
fprintf(fid,'%30.25f\n',y(end,:));
fclose(fid);

% MAXXO
% MAXY
 
mins = toc / 60;
mins

% --------------------------------------------------------------------------

function getamps(count)

global AMP IDENTAMPS NUMCOM

minmag = 0.25;
magstep = 0.5;
% Fig. 9D Edwards, Hirst & Suzuki
tt = [0.115 0.496 0.714 0.837 0.908 0.95 0.974 0.987 0.995];  % lammax = 28 / unscal
% Below it's extrapolated exponentially back towards zero (first bin is fuller)
% tt = [0.435 0.679 0.818 0.896 0.942 0.969 0.984 0.993 0.998];  % lammax = 43 / unscal

if (IDENTAMPS > 0)
    rand('state', 100000);    % include this call to get identical AMPs when
end                           %   PROSPECTS (or DURN) is varied

for ii = 1:NUMCOM
    AMP(ii,:) = rand(1,count);
    for kk = 1:count
        if (AMP(ii,kk) < tt(1))
            AMP(ii,kk) = minmag;
        elseif (AMP(ii,kk) < tt(2))
            AMP(ii,kk) = magstep + minmag;
        elseif (AMP(ii,kk) < tt(3))
            AMP(ii,kk) = 2 * magstep + minmag;
        elseif (AMP(ii,kk) < tt(4))
            AMP(ii,kk) = 3 * magstep + minmag;
        elseif (AMP(ii,kk) < tt(5))
            AMP(ii,kk) = 4 * magstep + minmag;
        elseif (AMP(ii,kk) < tt(6))
            AMP(ii,kk) = 5 * magstep + minmag;
        elseif (AMP(ii,kk) < tt(7))
            AMP(ii,kk) = 6 * magstep + minmag;
        elseif (AMP(ii,kk) < tt(8))
            AMP(ii,kk) = 7 * magstep + minmag;
        elseif (AMP(ii,kk) < tt(9))
            AMP(ii,kk) = 8 * magstep + minmag;
        else
            AMP(ii,kk) = 9 * magstep + minmag;
        end
    end
end

% --------------------------------------------------------------------------

function z = PlatCond(t,mess,comp)

global STARTS LAMBDA GAPS LAMINC NUMUNITS LGAPS LSTARTS AU BU KL
global GPLATM AMP RANDAMPS HILL PROSPECTS UNITSCAL MAXPROS PREV AMPINC

unitdur = 6.2 * UNITSCAL; % 3.1 same as 8.2; use 6.2; ICCim
                          % 8.2 works for diffexp GPLATM = 0.08, 0.453 & 0.083 sec.

ttemp = mess .^ HILL;
invlambda = (KL .^ HILL + ttemp) / (LAMINC * ttemp);
LAMBDA = 1.0 / invlambda;

for kk = 2:LGAPS(comp,:)              % could be Newton Raphson style (Quicksort)
    if (STARTS(comp,kk)>t)
        break
    end
end

if ((kk-PREV(comp)) > MAXPROS(comp))
    MAXPROS(comp) = kk - PREV(comp);
end
PREV(comp) = kk;

lastpros = kk + PROSPECTS;      % if (kk + PROSPECTS)>length(GAPS) maybe error
for i = kk:lastpros
    STARTS(comp,i) = STARTS(comp,i-1) + GAPS(comp,i) * invlambda;
end

early = t - unitdur;
for mm = 1:LSTARTS(comp,:)
    if (STARTS(comp,mm)>early)
        break
    end
end

z = 0;
for i = mm:LSTARTS(comp,:)
    t1 = t - STARTS(comp,i);
    if (t1>0)
        if (RANDAMPS > 0)
            z = z + AMP(comp,i) * (exp(AU * t1) - exp(BU * t1))^3;
        else
            z = z + (exp(AU * t1) - exp(BU * t1))^3;
        end
    else
        break
    end
end
z = GPLATM * (1 + AMPINC * LAMBDA / LAMINC) * z;

NUMUNITS = i-1;

% --------------------------------------------------------------------------

function dydt = f(t,y)

global PULON PULOFF IPULSE VM SM VH SH TAUM1 TAUH1 TAUHLIM TAUH2
global KF KB KM KH LAMINC STH HILL KL
global UNITS GPRIMAX TAUACT1 VACT SACT VINACT SINACT
global MBAS HBAS COMBAS COHBAS RICC TAUIN1 TAUIN2 TAUINLIM STIN
global ROUND NUMCOM IPULSEL PULONL PULOFFL PULONL2 PULOFFL2 NUDGE
global KFP KBP PEELED TRANS MAXXO MAXY

EPrim = -20; % -20;        % mV                              % -10 for decay overlap figure
EPlat = -20; % -20;        % -20 for CMB;          % mV      % -10 for decay overlap figure
ERest = -65;        % mV
gRest = 1.0 / RICC;  % nanoSiemens  % if RICC=12.3e-3, then gRest = 81.3nS;
Tm = 0.050;          % seconds      %0.019s in Goto, Matsuoka & Noma, 2004 (25.2pF 0.76gohm)
capac = Tm * gRest;   % nanoFarads

gLIC =  1000 / 3.27;  % 3.27;      %7.23 for small bundle;  % nanoSiemens
gL = 1000 / 9.28;     % 9.28;      %20.52 for small bundle; % nanoSiemens
TmL = 0.162;          % seconds     cousins et al 1993 say 140 ms in ilium
capacL = TmL * gL;    % nanoFarads  3-5 times larger gives better long muscle waveshape
EL = ERest;

if (ROUND > 0)
    gRest = 80;
    capac = 4;
    gLIC = 300;
    gL = 110;
    capacL = 20;
end

if (PEELED > 0)
    gLIC = 0;
end
gICIC = 35.0 * gRest;    %0.01  %35.0  %80.0   %%%%49     39     35     33     30
gLMLM = 52.4 * gL;       %60.0  %52.4  %0.01   %%%%49     52     52.4   53     54
                                               %%%%280ms  310ms  330ms  340ms  360ms (dodgy)
% gICIC = 0.1 * 35.0 * gRest;
% gLMLM = 0;

if (TRANS>0 & TRANS<1)
    gLMLM = TRANS * gLMLM;              % TRANS=0.2 ->       390ms  400ms
end
                                                        
for ii=1:NUMCOM

    voltindex = 9 * ii - 8;
    iccvolt = y(voltindex);
    
    if (ii > 1)
        vact1 = VACT;
        vinact1 = VINACT;
        tauinlim1 = TAUINLIM;
    else
        vact1 = VACT + NUDGE;
        vinact1 = VINACT + NUDGE;
        tauinlim1 = TAUINLIM + NUDGE;
    end
    
    gPrim = GPRIMAX * y(voltindex+7);

%  Primary component in 1st compartment only - careful
%     if (ii == 1)
%        gPrim = GPRIMAX * y(voltindex+7);
%     else
%        gPrim = 0;
%     end

    if (gPrim > MAXXO)
        MAXXO = gPrim;
        MAXY = y(voltindex+7);
    end
    
    actinf(ii) = 1 / (1 + exp((vact1 - iccvolt) * SACT));
    inactinf(ii) = 1 / (1 + exp((iccvolt - vinact1) * SINACT));
    
    tauact(ii) = TAUACT1;
    
    if (TAUIN1 ~= TAUIN2)
      tauinact(ii) = TAUIN1 + (TAUIN2 - TAUIN1) / (1 + exp((iccvolt - tauinlim1) * STIN));
    else
      tauinact(ii) = TAUIN1;
    end

    if (KF > 0)
        if (UNITS > 0)
            gPlat = PlatCond(t,y(voltindex+4),ii);
        else
            ttemp = y(voltindex+4) .^ HILL;
            gPlat = LAMINC * ttemp / (KL .^ HILL + ttemp);
        end
    else
        gPlat = 0;
    end
    
    minf(ii) = MBAS + COMBAS / (1 + exp((VM - iccvolt) * SM));
    hinf(ii) = HBAS + COHBAS / (1 + exp((iccvolt - VH) * SH));

    taumv(ii) = TAUM1;
      
    if (TAUH1 ~= TAUH2)
        tauhv(ii) = TAUH1 + (TAUH2 - TAUH1) / (1 + exp((iccvolt - TAUHLIM) * STH));
    else
        tauhv(ii) = TAUH1;
    end

    IPrim = gPrim * (iccvolt-EPrim);
    IPlat = gPlat * (iccvolt-EPlat);
    IRest = gRest * (iccvolt-ERest);

    lmvolt = y(voltindex+8);
        
    if (ii>1)
        IcicNear = IcicFar;
        IlmlmNear = IlmlmFar;
    else
        IcicNear = 0;
        IlmlmNear = 0;
    end
    if (ii<NUMCOM)
        IcicFar = gICIC * (iccvolt - y(voltindex+9));
        IlmlmFar = gLMLM * (lmvolt - y(voltindex+17));
    else
        IcicFar = 0;
        IlmlmFar = 0;
    end
 
    ILIC = gLIC * (iccvolt-lmvolt);
    Itot(ii) =  IPrim + IPlat + IRest + ILIC + IcicFar - IcicNear;
    Ilong(ii) = gL * (lmvolt-EL) - ILIC + IlmlmFar - IlmlmNear;
end

if (IPULSE~=0)
    if (t>PULON & t<PULOFF)
       Itot(1) = Itot(1) - IPULSE;
    end
end

if (IPULSEL~=0)
    if ((t>PULONL & t<PULOFFL) | (t>PULONL2 & t<PULOFFL2))
      Ilong(1) = Ilong(1) - IPULSEL;
    end
end

     dydt = [ -Itot(1) / capac
              (minf(1) - y(2)) / taumv(1)
              (hinf(1) - y(3)) / tauhv(1)
              KF * y(2) * y(3) - KB * y(4)
              KM * y(4) - KH * y(5)
              (actinf(1) - y(6)) / tauact(1)
              (inactinf(1) - y(7)) / tauinact(1)
              KFP * y(6) * y(7) - KBP * y(8)
              -Ilong(1) / capacL
              -Itot(2) / capac
              (minf(2) - y(11)) / taumv(2)
              (hinf(2) - y(12)) / tauhv(2)
              KF * y(11) * y(12) - KB * y(13)
              KM * y(13) - KH * y(14)
              (actinf(2) - y(15)) / tauact(2)
              (inactinf(2) - y(16)) / tauinact(2)
              KFP * y(15) * y(16) - KBP * y(17)
              -Ilong(2) / capacL
              -Itot(3) / capac
              (minf(3) - y(20)) / taumv(3)
              (hinf(3) - y(21)) / tauhv(3)
              KF * y(20) * y(21) - KB * y(22)
              KM * y(22) - KH * y(23)
              (actinf(3) - y(24)) / tauact(3)
              (inactinf(3) - y(25)) / tauinact(3)
              KFP * y(24) * y(25) - KBP * y(26)
              -Ilong(3) / capacL
              -Itot(4) / capac
              (minf(4) - y(29)) / taumv(4)
              (hinf(4) - y(30)) / tauhv(4)
              KF * y(29) * y(30) - KB * y(31)
              KM * y(31) - KH * y(32)
              (actinf(4) - y(33)) / tauact(4)
              (inactinf(4) - y(34)) / tauinact(4)
              KFP * y(33) * y(34) - KBP * y(35)
              -Ilong(4) / capacL
              -Itot(5) / capac
              (minf(5) - y(38)) / taumv(5)
              (hinf(5) - y(39)) / tauhv(5)
              KF * y(38) * y(39) - KB * y(40)
              KM * y(40) - KH * y(41)
              (actinf(5) - y(42)) / tauact(5)
              (inactinf(5) - y(43)) / tauinact(5)
              KFP * y(42) * y(43) - KBP * y(44)
              -Ilong(5) / capacL
              -Itot(6) / capac
              (minf(6) - y(47)) / taumv(6)
              (hinf(6) - y(48)) / tauhv(6)
              KF * y(47) * y(48) - KB * y(49)
              KM * y(49) - KH * y(50)
              (actinf(6) - y(51)) / tauact(6)
              (inactinf(6) - y(52)) / tauinact(6)
              KFP * y(51) * y(52) - KBP * y(53)
              -Ilong(6) / capacL
              -Itot(7) / capac
              (minf(7) - y(56)) / taumv(7)
              (hinf(7) - y(57)) / tauhv(7)
              KF * y(56) * y(57) - KB * y(58)
              KM * y(58) - KH * y(59)
              (actinf(7) - y(60)) / tauact(7)
              (inactinf(7) - y(61)) / tauinact(7)
              KFP * y(60) * y(61) - KBP * y(62)
              -Ilong(7) / capacL
              -Itot(8) / capac
              (minf(8) - y(65)) / taumv(8)
              (hinf(8) - y(66)) / tauhv(8)
              KF * y(65) * y(66) - KB * y(67)
              KM * y(67) - KH * y(68)
              (actinf(8) - y(69)) / tauact(8)
              (inactinf(8) - y(70)) / tauinact(8)
              KFP * y(69) * y(70) - KBP * y(71)
              -Ilong(8) / capacL
              -Itot(9) / capac
              (minf(9) - y(74)) / taumv(9)
              (hinf(9) - y(75)) / tauhv(9)
              KF * y(74) * y(75) - KB * y(76)
              KM * y(76) - KH * y(77)
              (actinf(9) - y(78)) / tauact(9)
              (inactinf(9) - y(79)) / tauinact(9)
              KFP * y(78) * y(79) - KBP * y(80)
              -Ilong(9) / capacL
              -Itot(10) / capac
              (minf(10) - y(83)) / taumv(10)
              (hinf(10) - y(84)) / tauhv(10)
              KF * y(83) * y(84) - KB * y(85)
              KM * y(85) - KH * y(86)
              (actinf(10) - y(87)) / tauact(10)
              (inactinf(10) - y(88)) / tauinact(10)
              KFP * y(87) * y(88) - KBP * y(89)
              -Ilong(10) / capacL
              -Itot(11) / capac
              (minf(11) - y(92)) / taumv(11)
              (hinf(11) - y(93)) / tauhv(11)
              KF * y(92) * y(93) - KB * y(94)
              KM * y(94) - KH * y(95)
              (actinf(11) - y(96)) / tauact(11)
              (inactinf(11) - y(97)) / tauinact(11)
              KFP * y(96) * y(97) - KBP * y(98)
              -Ilong(11) / capacL
              -Itot(12) / capac
              (minf(12) - y(101)) / taumv(12)
              (hinf(12) - y(102)) / tauhv(12)
              KF * y(101) * y(102) - KB * y(103)
              KM * y(103) - KH * y(104)
              (actinf(12) - y(105)) / tauact(12)
              (inactinf(12) - y(106)) / tauinact(12)
              KFP * y(105) * y(106) - KBP * y(107)
              -Ilong(12) / capacL
              -Itot(13) / capac
              (minf(13) - y(110)) / taumv(13)
              (hinf(13) - y(111)) / tauhv(13)
              KF * y(110) * y(111) - KB * y(112)
              KM * y(112) - KH * y(113)
              (actinf(13) - y(114)) / tauact(13)
              (inactinf(13) - y(115)) / tauinact(13)
              KFP * y(114) * y(115) - KBP * y(116)
              -Ilong(13) / capacL
              -Itot(14) / capac
              (minf(14) - y(119)) / taumv(14)
              (hinf(14) - y(120)) / tauhv(14)
              KF * y(119) * y(120) - KB * y(121)
              KM * y(121) - KH * y(122)
              (actinf(14) - y(123)) / tauact(14)
              (inactinf(14) - y(124)) / tauinact(14)
              KFP * y(123) * y(124) - KBP * y(125)
              -Ilong(14) / capacL
              -Itot(15) / capac
              (minf(15) - y(128)) / taumv(15)
              (hinf(15) - y(129)) / tauhv(15)
              KF * y(128) * y(129) - KB * y(130)
              KM * y(130) - KH * y(131)
              (actinf(15) - y(132)) / tauact(15)
              (inactinf(15) - y(133)) / tauinact(15)
              KFP * y(132) * y(133) - KBP * y(134)
              -Ilong(15) / capacL
              -Itot(16) / capac
              (minf(16) - y(137)) / taumv(16)
              (hinf(16) - y(138)) / tauhv(16)
              KF * y(137) * y(138) - KB * y(139)
              KM * y(139) - KH * y(140)
              (actinf(16) - y(141)) / tauact(16)
              (inactinf(16) - y(142)) / tauinact(16)
              KFP * y(141) * y(142) - KBP * y(143)
              -Ilong(16) / capacL
              -Itot(17) / capac
              (minf(17) - y(146)) / taumv(17)
              (hinf(17) - y(147)) / tauhv(17)
              KF * y(146) * y(147) - KB * y(148)
              KM * y(148) - KH * y(149)
              (actinf(17) - y(150)) / tauact(17)
              (inactinf(17) - y(151)) / tauinact(17)
              KFP * y(150) * y(151) - KBP * y(152)
              -Ilong(17) / capacL
              -Itot(18) / capac
              (minf(18) - y(155)) / taumv(18)
              (hinf(18) - y(156)) / tauhv(18)
              KF * y(155) * y(156) - KB * y(157)
              KM * y(157) - KH * y(158)
              (actinf(18) - y(159)) / tauact(18)
              (inactinf(18) - y(160)) / tauinact(18)
              KFP * y(159) * y(160) - KBP * y(161)
              -Ilong(18) / capacL
              -Itot(19) / capac
              (minf(19) - y(164)) / taumv(19)
              (hinf(19) - y(165)) / tauhv(19)
              KF * y(164) * y(165) - KB * y(166)
              KM * y(166) - KH * y(167)
              (actinf(19) - y(168)) / tauact(19)
              (inactinf(19) - y(169)) / tauinact(19)
              KFP * y(168) * y(169) - KBP * y(170)
              -Ilong(19) / capacL
              -Itot(20) / capac
              (minf(20) - y(173)) / taumv(20)
              (hinf(20) - y(174)) / tauhv(20)
              KF * y(173) * y(174) - KB * y(175)
              KM * y(175) - KH * y(176)
              (actinf(20) - y(177)) / tauact(20)
              (inactinf(20) - y(178)) / tauinact(20)
              KFP * y(177) * y(178) - KBP * y(179)
              -Ilong(20) / capacL
            ];

% --------------------------------------------------------------------------

function inacts
     
global VACT SACT VINACT SINACT

i = 1:71;
em(i) = i - 91;
act(i) = 1 ./ (1 + exp((VACT - em) * SACT));
inact(i) = 1 ./ (1 + exp((em - VINACT) * SINACT));
mexp(i) = act(i) .*act(i);
fac(i) = mexp(i) .* inact(i);

figure(202);
plot(em,mexp,em,inact,em,fac,em,act);
axis([-65 -45 0 1]);
% --------------------------------------------------------------------------
