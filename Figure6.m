function batcher
%
% 14mm/s circumferential propagation; 
% 600um CM compartment;
% 43ms CM compartment transit time;
% 3mm is CM length constant
%

clear global;

global DURN STARTS LAMBDA GAPS PULON PULOFF IPULSE VM SM VH SH
global LAMINC TAUH1 TAUHLIM TAUH2 NUMUNITS UNITS DETAILS NEWNOISE
global KF KB KM KH KL STH FIGNUM GPRIMAX TAUACT1 VACT SACT VINACT SINACT
global MBAS HBAS VTOP VBOT AU BU GPLATM AMP RANDAMPS HILL RICC AX1
global COMBAS COHBAS TBOT LGAPS LSTARTS
global TAUIN1 STIN PROSPECTS MAXINTSTEP RATEFAC UNITSCAL
global MAXPROS PREV TAUM1 STM IDENTAMPS
global PROSPECTSC LAMINCC LAMBDAC NUMUNITSC GAPSC LGAPSC STARTSC LSTARTSC
global AMPC KFC SMC SHC MAXPROSC PREVC
global AMPINC AMPINCC
global PEELED TAUH1C NUMCOM PASSIVE TRANS
global ROUND TAUACT1 TAUINLIM TAUIN2 STANDICS LBOT LTOP KFP KBP
global NUMICCMY NUMICCUSED NUDGE TAUH2C
global GAUSSTCL HY0 HA HOM HXC TAUNUM

tic;

NUMCOM = 16;
NUMICCMY = 6;
NUMICCUSED = 6; %2; %6;
PASSIVE = 1; %1; %0;
PEELED = 0;
TRANS = 0.1;        %0.1
ROUND = 1;
DETAILS = 0;
DURN = 51; %15; %70; %60; %34;      % 1200 works ok on Dell 3GHz;
FIGNUM = 20;
TAUNUM = 200;
MAXINTSTEP = 3e-3; %6e-3; %3e-3; %6e-4 differs from 9e-4;  % sec
PROSPECTS = 20;   %18;
IDENTAMPS = 0;          % can be set when PROSPECTS or DURN is varied
STANDICS = 0;           % -1 reads/writes FCs; 0 reads FCs; 1 uses equal values; 2 uses equal values and writes FCs;
VTOP = 0;
VBOT = -80;
TBOT = 0;
LBOT = 2.5;
LTOP = 4.0;
%
%   Unit parameters
%
UNITS = 1; %0; %1;      % if UNITS = 0, use RATEFAC = 12?
RANDAMPS = 1;
NEWNOISE = 0; %0;
RICC = 12.3e-3; %27.2e-3 for small bundle;      % gigohms  see Cousins et al 2003, Table 1
%       12.3e-3 gives ~100uV minimum unit in CMB, ~300uV in ICCmy
UNITSCAL = 1; %1;       % multiplies unit size, RICC independent
RATEFAC = 2; %2; %12; %2;    % 2.7 LONGYES = 0;  % 1 for CMB; multiplies max frequency
if (UNITS == 0)
    RATEFAC = 6 * RATEFAC;
end

unitbase = 0.2; %0.2;          % 0.2,25 gives 0.5mV unit
unscal = UNITSCAL * unitbase;
GPLATM = unscal / RICC;  % multiplies unit size (see  above)
AU = -1.0 / 0.434;       % average of 17 Edw, Hir & Suz
BU = -1.0 / 0.077;       % average of 17 Edw, Hir & Suz
%
%   primary component voltage dependencies & kinetics
%
vshift = 0;             %0; % -6 for spon
NUDGE = -10; %-10; %-10.7; %-11.2;             %-8; %-9 slows velocity & speeds rate
GPRIMAX = 1.2*2.035e7;      %2.035e7   <<<<<<<<<<<<<<<<

KFP = 0.012;            %0.012;    %0.0114 slows velocity & rate;
KBP = 6;                %6;        %6.5 speeds velocity and slows rate

TAUACT1 = 0.06;             %0.06;     % 0.06 gives pointier ears than 0.1; ok for UNITS=1
TAUIN1 = 0.25;              %0.25;
TAUINLIM = -48 + vshift;    %-48 + vshift;
STIN = 2;                   %2;
TAUIN2 = 90;                %90;            % 3 for Hiroe's spikes during diastole
VACT = -48 + vshift;        %-48 + vshift;
SACT = 0.78;                %0.78;           % 0.9 slows rate, 0.7 speeds rate
VINACT = -53 + vshift;      %-53 + vshift;   % -50 speeds rate, -56 slows rate 
SINACT = 1;                 %1;              % 2 speeds rate
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

TAUM1 = 0.15;               %0.15;
TAUH1 = 1.7;                %1.7; %2.1; %2.5; %3;   %1.2 in CMB paper
TAUHLIM = -45.5;            %-45.5
STH = 0.5;                  %0.5
TAUH2 = 6;                  %6

GAUSSTCL = 0;
HY0 = TAUH1;
HA = TAUH2 - TAUH1;
HOM = 6;        %12;
HXC = VH+15;    %VH;

%   current injection into CMB
IPULSE = 0; %1e3;
PULON = 9999; %38;
PULOFF = 10000; %300;

COMBAS = 1 - MBAS;
COHBAS = 1 - HBAS;
NUMUNITS = 0;
MAXPROS = 0;
PREV = 0;

    PROSPECTSC = 12;                % 12;
    LAMINCC = 1.0 * LAMINC;         % 1.0 * LAMINC;  % 15*28/unscal;
    LAMBDAC = 1;                    % 1;
    AMPINCC = (AMPINC+1) / 2 - 1;
    KFC = 2 * 1.25 * 0.0157;      % 1.3 * 1.25 * 0.0157;  <<<<<<<<<
    SMC = 0.3;                      % 0.3;
    SHC = 0.32;                     % 0.32;
    TAUH1C = TAUH1;
    TAUH2C = TAUH2;
    NUMUNITSC = 0;
    for ii = 1:NUMCOM
        MAXPROSC(ii) = 0;
        PREVC(ii) = 0;
    end

if (PASSIVE > 0)
    PEELED = 1;
    IPULSE = 100000;
    PULON = 5;
    PULOFF = 10;
    DURN = 10;
    KFC = 0;
    STANDICS = 1;           % 1 uses equal values;
end

if (NEWNOISE == 0)
   rand('state', 0);
end

d1i;

if (NEWNOISE ~= 0)
    NEWNOISE
end

if (GAUSSTCL > 0)
    GAUSSTCL
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
global DURN STARTS LAMBDA GAPS IPULSE LAMINC NUMUNITS KL FIGNUM
global VTOP VBOT TBOT LGAPS LSTARTS NUMCOM UNITS DETAILS PASSIVE
global RANDAMPS HILL PROSPECTS MAXINTSTEP RATEFAC MAXPROS MAXPROSC
global PROSPECTSC LAMINCC LAMBDAC NUMUNITSC GAPSC LGAPSC STARTSC LSTARTSC
global STANDICS LBOT LTOP

vi0 = -65;         % -65;
m0 = 0.04;         % 0.04;
h0 = 0.001;        % 0.001;
interm0 = 40e-6;   % 40e-6;
mess0 = 200e-6;    % 200e-6;
act0 = 0.01;       % 0.01;        % 0.0033 for CMB;
inact0 = 0.001;    % 0.001;       % 0.2489 for CMB;
messP0 = 1e-7;     % 1e-7;
vLong0 = -65;      % -65;
mC0 = 0.18;        %0.18;
hC0 = 0.028;       %0.28;
intermC0 = 210e-6; %210e-6;
messC0 = 340e-6;   %340e-6;
vCMB0 = -65;       % -65;
actC0 = 0.01;      % 0.01;
inactC0 = 0.001;   % 0.001;

LAMBDA = 0.3;      % 0.3;
evlen = 0;
evlenC = 0;

% plattau;

if (DETAILS == 1)
    figure(FIGNUM);
end

    if (STANDICS > 0)
        y0 = [mC0;hC0;intermC0;messC0;vCMB0;mC0;hC0;intermC0;messC0;vCMB0;...
              mC0;hC0;intermC0;messC0;vCMB0;mC0;hC0;intermC0;messC0;vCMB0;...
              mC0;hC0;intermC0;messC0;vCMB0;mC0;hC0;intermC0;messC0;vCMB0;...
              mC0;hC0;intermC0;messC0;vCMB0;mC0;hC0;intermC0;messC0;vCMB0;...
              mC0;hC0;intermC0;messC0;vCMB0;mC0;hC0;intermC0;messC0;vCMB0;...
              mC0;hC0;intermC0;messC0;vCMB0;mC0;hC0;intermC0;messC0;vCMB0;...
              mC0;hC0;intermC0;messC0;vCMB0;mC0;hC0;intermC0;messC0;vCMB0;...
              mC0;hC0;intermC0;messC0;vCMB0;mC0;hC0;intermC0;messC0;vCMB0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0;...
              vi0;m0;h0;interm0;mess0;act0;inact0;messP0;vLong0...
            ];
    else
%        fid = fopen('C:\Matfiles\FCs134.dat');
        fid = fopen('FCs134.dat');
        y0 = fscanf(fid,'%g',[134 inf]);    % It has 134 rows.
        y0 = y0';
        fclose(fid);
    end

if (UNITS > 0)
    if (DURN <= 20)
        evlen = round(LAMINC * DURN + PROSPECTS);
    else
        evlen = round(0.6 * LAMINC * DURN + PROSPECTS);   % 0.6
    end
    GAPS = -log(rand(1,evlen));
    LGAPS = length(GAPS);
    STARTS = cumsum(GAPS / LAMBDA);
    LSTARTS = length(STARTS);

        if (DURN <= 20)
            evlenC = round(LAMINCC * DURN + PROSPECTSC);
        else
            evlenC = round(0.4 * LAMINCC * DURN + PROSPECTSC); %0.4
        end

    for ii = 1:NUMCOM
        GAPSC(ii,:) = -log(rand(1,evlenC));
        LGAPSC(ii,:) = length(GAPSC(ii,:));
        STARTSC(ii,:) = cumsum(GAPSC(ii,:) / LAMBDAC);
        LSTARTSC(ii,:) = length(STARTSC(ii,:));
    end

    if (RANDAMPS > 0)
        getamps(evlen,evlenC);
    end
end

    options = [];
    options = odeset('MaxStep', MAXINTSTEP);
    [t,y] = ode15s(@f,[0 DURN],y0,options);

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
% this plots activation & inactivation for Primary component
    figure(270);
    plot(t,y(:,6),t,y(:,7),t,y(:,6) .* y(:,7));
   axis([TBOT DURN -0.05 1.05]);
% circular muscle
    figure(271);
    plot(t,y(:,9),t,y(:,10));
    axis([TBOT DURN 0 1]);
    figure(272);
    plot(t,y(:,11),t,y(:,12));
    axis([TBOT DURN 0 0.005]);
    ttemp = y(:,12) .^ HILL;
    ttemp2 = KL ^ HILL;
    urate = (LAMINCC * ttemp) ./ (ttemp2 + ttemp);
    figure(273);
    plot(t,urate);
    axis([TBOT DURN 0 150*RATEFAC]);
end

figure(FIGNUM);
    plot(t,y(:,5),t,y(:,10),t,y(:,15),t,y(:,20),t,y(:,25),t,y(:,30),t,y(:,35),t,y(:,40),t,y(:,45),t,y(:,50),t,y(:,55),t,y(:,60),t,y(:,65),t,y(:,70),t,y(:,75),t,y(:,80),...
        t,y(:,81),t,y(:,90),t,y(:,99),t,y(:,108),t,y(:,117),t,y(:,126));
%    plot(t,y(:,81),t,y(:,90),t,y(:,99),t,y(:,108),t,y(:,117),t,y(:,126));
    axis([TBOT DURN VBOT VTOP]);
figure(FIGNUM+1);
%     plot(t,y(:,5),t,y(:,10),t,y(:,15),t,y(:,20),t,y(:,25),t,y(:,30),t,y(:,35),t,y(:,40),t,y(:,45),t,y(:,50),t,y(:,55),t,y(:,60),t,y(:,65),t,y(:,70),t,y(:,75),t,y(:,80));
if (PASSIVE > 0)
    axis([TBOT DURN VBOT VTOP]);
    len = length(y(:,1));
    prepwidth = 0.1;
    preplength = 0.6;
    if (IPULSE~=0)
        for i = 1:NUMCOM-1
            over = 5 * (i + 1);
            under = over - 5;
            lcCM(i) = preplength / (-log((y(len-1,over) - y(2,over)) / (y(len-1,under) - y(2,under))));
            lcCMone(i) = preplength * i / (-log((y(len-1,over) - y(2,over)) / (y(len-1,5) - y(2,5))));
            pos(i) = i;
        end
%        figure(FIGNUM+2);
        plot(pos,lcCM,pos,lcCMone);
        grid on
        axis([0 NUMCOM LBOT LTOP]);
    end
else
%    axis([DURN-20 DURN -48 -30]);
%    figure(FIGNUM+2);
    plot(t,y(:,5),t,y(:,55));
%    axis([DURN-20 DURN -48 -30]);
    axis([TBOT DURN VBOT VTOP]);
end

if (NUMUNITS > evlen)
    NUMUNITS
    evlen
end
if (NUMUNITSC > evlenC)
    NUMUNITSC
    evlenC
end
if (MAXPROS > PROSPECTS)
    MAXPROS
    PROSPECTS
end    
for ii = 1:NUMCOM
    if (MAXPROSC(ii) > PROSPECTSC)
        ii
        MAXPROSC(ii)
        PROSPECTSC
    end    
end

if ((STANDICS < 0) | (STANDICS == 2))
%    fid = fopen('C:\Matfiles\FCs134.dat','w');
    fid = fopen('FCs134.dat','w');
    fprintf(fid,'%30.25f\n',y(end,:));
    fclose(fid);
end

mins = toc / 60;
mins

% --------------------------------------------------------------------------

function getamps(count,countC)

global AMP IDENTAMPS AMPC NUMCOM

minmag = 0.25;
magstep = 0.5;
% Fig. 9D Edwards, Hirst & Suzuki
tt = [0.115 0.496 0.714 0.837 0.908 0.95 0.974 0.987 0.995];  % lammax = 28 / unscal
% Below it's extrapolated exponentially back towards zero (first bin is fuller)
% tt = [0.435 0.679 0.818 0.896 0.942 0.969 0.984 0.993 0.998];  % lammax = 43 / unscal

if (IDENTAMPS > 0)
    rand('state', 100000);    % include this call to get identical AMPs when
end                           %   PROSPECTS (or DURN) is varied

AMP = rand(1,count);

for kk = 1:count
    if (AMP(kk) < tt(1))
        AMP(kk) = minmag;
    elseif (AMP(kk) < tt(2))
        AMP(kk) = magstep + minmag;
    elseif (AMP(kk) < tt(3))
        AMP(kk) = 2 * magstep + minmag;
    elseif (AMP(kk) < tt(4))
        AMP(kk) = 3 * magstep + minmag;
    elseif (AMP(kk) < tt(5))
        AMP(kk) = 4 * magstep + minmag;
    elseif (AMP(kk) < tt(6))
        AMP(kk) = 5 * magstep + minmag;
    elseif (AMP(kk) < tt(7))
        AMP(kk) = 6 * magstep + minmag;
    elseif (AMP(kk) < tt(8))
        AMP(kk) = 7 * magstep + minmag;
    elseif (AMP(kk) < tt(9))
        AMP(kk) = 8 * magstep + minmag;
    else
        AMP(kk) = 9 * magstep + minmag;
    end
end

for ii = 1:NUMCOM
    AMPC(ii,:) = rand(1,countC);
    for kk = 1:countC
        if (AMPC(ii,kk) < tt(1))
            AMPC(ii,kk) = minmag;
        elseif (AMPC(ii,kk) < tt(2))
            AMPC(ii,kk) = magstep + minmag;
        elseif (AMPC(ii,kk) < tt(3))
            AMPC(ii,kk) = 2 * magstep + minmag;
        elseif (AMPC(ii,kk) < tt(4))
            AMPC(ii,kk) = 3 * magstep + minmag;
        elseif (AMPC(ii,kk) < tt(5))
            AMPC(ii,kk) = 4 * magstep + minmag;
        elseif (AMPC(ii,kk) < tt(6))
            AMPC(ii,kk) = 5 * magstep + minmag;
        elseif (AMPC(ii,kk) < tt(7))
            AMPC(ii,kk) = 6 * magstep + minmag;
        elseif (AMPC(ii,kk) < tt(8))
            AMPC(ii,kk) = 7 * magstep + minmag;
        elseif (AMPC(ii,kk) < tt(9))
            AMPC(ii,kk) = 8 * magstep + minmag;
        else
            AMPC(ii,kk) = 9 * magstep + minmag;
        end
    end
end
% --------------------------------------------------------------------------

function z = PlatCond(t,mess)

global STARTS LAMBDA GAPS LAMINC NUMUNITS LGAPS LSTARTS AU BU
global GPLATM AMP RANDAMPS KL HILL PROSPECTS UNITSCAL MAXPROS PREV AMPINC

unitdur = 6.2 * UNITSCAL; % 3.1 same as 8.2; use 6.2; ICCmy
                          % 8.2 works for diffexp GPLATM = 0.08, 0.453 & 0.083 sec.

ttemp = mess .^ HILL;
invlambda = (KL .^ HILL + ttemp) / (LAMINC * ttemp);
LAMBDA = 1.0 / invlambda;

for kk = 2:LGAPS              % could be Newton Raphson style (Quicksort)
    if (STARTS(kk)>t)
        break
    end
end

if ((kk-PREV) > MAXPROS)
    MAXPROS = kk - PREV;
end
PREV = kk;

lastpros = kk + PROSPECTS;      % if (kk + PROSPECTS)>length(GAPS) maybe error
for i = kk:lastpros
    STARTS(i) = STARTS(i-1) + GAPS(i) * invlambda;
end

early = t - unitdur;
for mm = 1:LSTARTS
    if (STARTS(mm)>early)
        break
    end
end

z = 0;
for i = mm:LSTARTS
    t1 = t - STARTS(i);
    if (t1>0)
        if (RANDAMPS > 0)
            z = z + AMP(i) * (exp(AU * t1) - exp(BU * t1))^3;
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

function z = CMBCond(t,mess,comp)

global STARTSC LAMBDAC GAPSC LAMINCC NUMUNITSC LGAPSC LSTARTSC AU BU
global GPLATM AMPC RANDAMPS KL HILL PROSPECTSC UNITSCAL MAXPROSC PREVC AMPINCC

unitdur = 6.2 * UNITSCAL; % 3.1 same as 8.2; use 6.2; ICCim
                          % 8.2 works for diffexp GPLATM = 0.08, 0.453 & 0.083 sec.

ttemp = mess .^ HILL;
invlambda = (KL .^ HILL + ttemp) / (LAMINCC * ttemp);
LAMBDAC = 1.0 / invlambda;

for kk = 2:LGAPSC(comp,:)              % could be Newton Raphson style (Quicksort)
    if (STARTSC(comp,kk)>t)
        break
    end
end

if ((kk-PREVC(comp)) > MAXPROSC(comp))
    MAXPROSC(comp) = kk - PREVC(comp);
end
PREVC(comp) = kk;

lastpros = kk + PROSPECTSC;      % if (kk + PROSPECTS)>length(GAPS) maybe error
for i = kk:lastpros
    STARTSC(comp,i) = STARTSC(comp,i-1) + GAPSC(comp,i) * invlambda;
end

early = t - unitdur;
for mm = 1:LSTARTSC(comp,:)
    if (STARTSC(comp,mm)>early)
        break
    end
end

z = 0;
for i = mm:LSTARTSC(comp,:)
    t1 = t - STARTSC(comp,i);
    if (t1>0)
        if (RANDAMPS > 0)
            z = z + AMPC(comp,i) * (exp(AU * t1) - exp(BU * t1))^3;
        else
            z = z + (exp(AU * t1) - exp(BU * t1))^3;
        end
    else
        break
    end
end
z = GPLATM * (1 + AMPINCC * LAMBDAC / LAMINCC) * z;

NUMUNITSC = i-1;

% --------------------------------------------------------------------------

function dydt = f(t,y)

global PULON PULOFF IPULSE VM SM VH SH TAUH1 TAUHLIM TAUH2
global KF KB KM KH LAMINC STH HILL KL
global UNITS GPRIMAX TAUACT1 VACT SACT VINACT SINACT
global MBAS HBAS COMBAS COHBAS RICC TAUIN1 TAUIN2 TAUINLIM STIN 
global TAUM1 KFC SMC SHC LAMINCC KFP KBP
global ROUND PEELED TAUH1C NUMCOM MESSOFF VOLTOFF NEXTOFF NUMICCMY TRANS
global NUMICCUSED NUDGE TAUH2C
global GAUSSTCL HY0 HA HOM HXC


EPrim = -20; % -20;        % mV                              % -10 for decay overlap figure
EPlat = -20; % -20;        % -20 for CMB;          % mV      % -10 for decay overlap figure
ERest = -65;        % mV
gRest = 1.0 / RICC;  % nanoSiemens
Tm = 0.050;          % seconds      %0.019s in Goto, Matsuoka & Noma, 2004 (25.2pF 0.76gohm)
capac = Tm * gRest;   % nanoFarads

gLIC =  1000 / 3.27;  % 3.27;      %7.23 for small bundle;  % nanoSiemens
gL = 1000 / 9.28;     % 9.28;      %20.52 for small bundle; % nanoSiemens
TmL = 0.162;          % seconds     cousins et al 1993 say 140 ms in ilium
capacL = TmL * gL;    % nanoFarads  3-5 times larger gives better long muscle waveshape
EL = ERest;

gC = 1000 / 1.99;
TmC = 0.162;         % seconds      also Table I
capacC = TmC * gC;   % nanoFarads
EC = ERest;

if (ROUND > 0)
    gRest = 80;
    capac = 4;
    gLIC = 300;
    gL = 110;
    capacL = 20;
    gC = 500;
    capacC = 80;
end

gCIC = gLIC;        % 0.1*gLIC, 10*gLIC, 100*gLIC, 1000*gLIC
if (PEELED > 0)
    gCIC = 1e-20 * gCIC;
end

gCC = 23.8 * gC;            % 23.8*gC gives 3mm length const
gICIC = 35.0 * gRest;       %0.01  %35.0  %80.0   %%%%49     39     35     33     30
gLMLM = TRANS * 52.4 * gL;  %60.0  %52.4  %0.01   %%%%49     52     52.4   53     54

% gCIC = 0;
% gLIC = 0;
% gCC = 0;
% gICIC = 0;
% gLMLM = 0;

for ii=1:NUMCOM

    voltindex = 5 * ii;
    cmvolt = y(voltindex);
    
    if (KFC > 0)
        if (UNITS > 0)
            gCMB(ii) = CMBCond(t,y(voltindex-1),ii);
        else
            ttemp = y(voltindex-1) .^ HILL;
            gCMB(ii) = LAMINCC * ttemp / (KL .^ HILL + ttemp);
        end
    else
        gCMB(ii) = 0;
    end
    
    minfC(ii) = MBAS + COMBAS / (1 + exp((VM - cmvolt) * SMC));
    hinfC(ii) = HBAS + COHBAS / (1 + exp((cmvolt - VH) * SHC));

    taumvC(ii) = TAUM1;
       
%     if (TAUH1C ~= TAUH2C)
%         tauhvC(ii) = TAUH1C + (TAUH2C - TAUH1C) / (1 + exp((cmvolt - TAUHLIM) * STH));
%     else
%         tauhvC(ii) = TAUH1C;
%     end

    if (GAUSSTCL > 0)
        tauhvC(ii) = HY0 + HA  .* exp(-0.5 .* ((cmvolt-HXC)./HOM).^2);
    else
        if (TAUH1C ~= TAUH2C)
            tauhvC(ii) = TAUH1C + (TAUH2C - TAUH1C) / (1 + exp((cmvolt - TAUHLIM) * STH));
        else
            tauhvC(ii) = TAUH1C;
        end
    end
    
    if (ii>1)
        ICmbCmbNear = ICmbCmbFar;
    else
        ICmbCmbNear = 0;
    end
    if (ii<NUMCOM)
        ICmbCmbFar = gCC * (cmvolt - y(voltindex+5));
    else
        ICmbCmbFar = 0;
    end

    Icmb(ii) = gC*(cmvolt-EC) + gCMB(ii)*(cmvolt-EPlat) + ICmbCmbFar - ICmbCmbNear;

end

for ii=1:NUMICCMY

    sindex = 81 + 9 * (ii-1);
    iccvolt = y(sindex);

    if (ii > 1)
        vact1 = VACT;
        vinact1 = VINACT;
        tauinlim1 = TAUINLIM;
    else
        vact1 = VACT + NUDGE;
        vinact1 = VINACT + NUDGE;
        tauinlim1 = TAUINLIM + NUDGE;
    end

    gPrim = GPRIMAX * y(sindex+7);

    actinf(ii) = 1 / (1 + exp((vact1 - iccvolt) * SACT));
    inactinf(ii) = 1 / (1 + exp((iccvolt - vinact1) * SINACT));
    
    tauact(ii) = TAUACT1;
    
    if (TAUIN1 ~= TAUIN2)
      tauinact(ii) = TAUIN1 + (TAUIN2 - TAUIN1) / (1 + exp((iccvolt - tauinlim1) * STIN));
    else
      tauinact(ii) = TAUIN1;
    end

    if (UNITS > 0)
        gPlat = PlatCond(t,y(sindex+4));
    else
        ttemp = y(sindex+4) .^ HILL;
        gPlat = LAMINC * ttemp / (KL .^ HILL + ttemp);
    end
    
    minf(ii) = MBAS + COMBAS / (1 + exp((VM - iccvolt) * SM));
    hinf(ii) = HBAS + COHBAS / (1 + exp((iccvolt - VH) * SH));

    taumv(ii) = TAUM1;
      
%     if (TAUH1 ~= TAUH2)
%         tauhv(ii) = TAUH1 + (TAUH2 - TAUH1) / (1 + exp((iccvolt - TAUHLIM) * STH));
%     else
%         tauhv(ii) = TAUH1;
%     end

   if (GAUSSTCL > 0)
        tauhv(ii) = HY0 + HA  .* exp(-0.5 .* ((iccvolt-HXC)./HOM).^2);
   else
       if (TAUH1 ~= TAUH2)
           tauhv(ii) = TAUH1 + (TAUH2 - TAUH1) / (1 + exp((iccvolt - TAUHLIM) * STH));
       else
           tauhv(ii) = TAUH1;
       end
   end

    IPrim = gPrim * (iccvolt-EPrim);
    IPlat = gPlat * (iccvolt-EPlat);
    IRest = gRest * (iccvolt-ERest);

    lmvolt = y(sindex+8);
        
    if (ii>1)
        IcicNear = IcicFar;
        IlmlmNear = IlmlmFar;
    else
        IcicNear = 0;
        IlmlmNear = 0;
    end
    if (ii<NUMICCMY)
        IcicFar = gICIC * (iccvolt - y(sindex+9));
        IlmlmFar = gLMLM * (lmvolt - y(sindex+17));
    else
        IcicFar = 0;
        IlmlmFar = 0;
    end
 
    ILIC = gLIC * (iccvolt-lmvolt);
    Ilong(ii) = gL * (lmvolt-EL) - ILIC + IlmlmFar - IlmlmNear;
    Itot(ii) = IPrim + IPlat + IRest + ILIC + IcicFar - IcicNear;

    if (ii<=NUMICCUSED)
        ICIC = gCIC * (iccvolt-y(5*ii));
        Itot(ii) = Itot(ii) + ICIC;
        Icmb(ii) = Icmb(ii) - ICIC;
    end

end

if (IPULSE~=0)
    if (t>PULON & t<PULOFF)
       Icmb(1) = Icmb(1) - IPULSE;
    end
end

     dydt = [ (minfC(1) - y(1)) / taumvC(1)
              (hinfC(1) - y(2)) / tauhvC(1)
              KFC * y(1) * y(2) - KB * y(3)
              KM * y(3) - KH * y(4)
              -Icmb(1) / capacC
              (minfC(2) - y(6)) / taumvC(2)
              (hinfC(2) - y(7)) / tauhvC(2)
              KFC * y(6) * y(7) - KB * y(8)
              KM * y(8) - KH * y(9)
              -Icmb(2) / capacC
              (minfC(3) - y(11)) / taumvC(3)
              (hinfC(3) - y(12)) / tauhvC(3)
              KFC * y(11) * y(12) - KB * y(13)
              KM * y(13) - KH * y(14)
              -Icmb(3) / capacC
              (minfC(4) - y(16)) / taumvC(4)
              (hinfC(4) - y(17)) / tauhvC(4)
              KFC * y(16) * y(17) - KB * y(18)
              KM * y(18) - KH * y(19)
              -Icmb(4) / capacC
              (minfC(5) - y(21)) / taumvC(5)
              (hinfC(5) - y(22)) / tauhvC(5)
              KFC * y(21) * y(22) - KB * y(23)
              KM * y(23) - KH * y(24)
              -Icmb(5) / capacC
              (minfC(6) - y(26)) / taumvC(6)
              (hinfC(6) - y(27)) / tauhvC(6)
              KFC * y(26) * y(27) - KB * y(28)
              KM * y(28) - KH * y(29)
              -Icmb(6) / capacC
              (minfC(7) - y(31)) / taumvC(7)
              (hinfC(7) - y(32)) / tauhvC(7)
              KFC * y(31) * y(32) - KB * y(33)
              KM * y(33) - KH * y(34)
              -Icmb(7) / capacC
              (minfC(8) - y(36)) / taumvC(8)
              (hinfC(8) - y(37)) / tauhvC(8)
              KFC * y(36) * y(37) - KB * y(38)
              KM * y(38) - KH * y(39)
              -Icmb(8) / capacC
              (minfC(9) - y(41)) / taumvC(9)
              (hinfC(9) - y(42)) / tauhvC(9)
              KFC * y(41) * y(42) - KB * y(43)
              KM * y(43) - KH * y(44)
              -Icmb(9) / capacC
              (minfC(10) - y(46)) / taumvC(10)
              (hinfC(10) - y(47)) / tauhvC(10)
              KFC * y(46) * y(47) - KB * y(48)
              KM * y(48) - KH * y(49)
              -Icmb(10) / capacC
              (minfC(11) - y(51)) / taumvC(11)
              (hinfC(11) - y(52)) / tauhvC(11)
              KFC * y(51) * y(52) - KB * y(53)
              KM * y(53) - KH * y(54)
              -Icmb(11) / capacC
              (minfC(12) - y(56)) / taumvC(12)
              (hinfC(12) - y(57)) / tauhvC(12)
              KFC * y(56) * y(57) - KB * y(58)
              KM * y(58) - KH * y(59)
              -Icmb(12) / capacC
              (minfC(13) - y(61)) / taumvC(13)
              (hinfC(13) - y(62)) / tauhvC(13)
              KFC * y(61) * y(62) - KB * y(63)
              KM * y(63) - KH * y(64)
              -Icmb(13) / capacC
              (minfC(14) - y(66)) / taumvC(14)
              (hinfC(14) - y(67)) / tauhvC(14)
              KFC * y(66) * y(67) - KB * y(68)
              KM * y(68) - KH * y(69)
              -Icmb(14) / capacC
              (minfC(15) - y(71)) / taumvC(15)
              (hinfC(15) - y(72)) / tauhvC(15)
              KFC * y(71) * y(72) - KB * y(73)
              KM * y(73) - KH * y(74)
              -Icmb(15) / capacC
              (minfC(16) - y(76)) / taumvC(16)
              (hinfC(16) - y(77)) / tauhvC(16)
              KFC * y(76) * y(77) - KB * y(78)
              KM * y(78) - KH * y(79)
              -Icmb(16) / capacC
              -Itot(1) / capac
              (minf(1) - y(82)) / taumv(1)
              (hinf(1) - y(83)) / tauhv(1)
              KF * y(82) * y(83) - KB * y(84)
              KM * y(84) - KH * y(85)
              (actinf(1) - y(86)) / tauact(1)
              (inactinf(1) - y(87)) / tauinact(1)
              KFP * y(86) * y(87) - KBP * y(88)
              -Ilong(1) / capacL
              -Itot(2) / capac
              (minf(2) - y(91)) / taumv(2)
              (hinf(2) - y(92)) / tauhv(2)
              KF * y(91) * y(92) - KB * y(93)
              KM * y(93) - KH * y(94)
              (actinf(2) - y(95)) / tauact(2)
              (inactinf(2) - y(96)) / tauinact(2)
              KFP * y(95) * y(96) - KBP * y(97)
              -Ilong(2) / capacL
              -Itot(3) / capac
              (minf(3) - y(100)) / taumv(3)
              (hinf(3) - y(101)) / tauhv(3)
              KF * y(100) * y(101) - KB * y(102)
              KM * y(102) - KH * y(103)
              (actinf(3) - y(104)) / tauact(3)
              (inactinf(3) - y(105)) / tauinact(3)
              KFP * y(104) * y(105) - KBP * y(106)
              -Ilong(3) / capacL
              -Itot(4) / capac
              (minf(4) - y(109)) / taumv(4)
              (hinf(4) - y(110)) / tauhv(4)
              KF * y(109) * y(110) - KB * y(111)
              KM * y(111) - KH * y(112)
              (actinf(4) - y(113)) / tauact(4)
              (inactinf(4) - y(114)) / tauinact(4)
              KFP * y(113) * y(114) - KBP * y(115)
              -Ilong(4) / capacL
              -Itot(5) / capac
              (minf(5) - y(118)) / taumv(5)
              (hinf(5) - y(119)) / tauhv(5)
              KF * y(118) * y(119) - KB * y(120)
              KM * y(120) - KH * y(121)
              (actinf(5) - y(122)) / tauact(5)
              (inactinf(5) - y(123)) / tauinact(5)
              KFP * y(122) * y(123) - KBP * y(124)
              -Ilong(5) / capacL
              -Itot(6) / capac
              (minf(6) - y(127)) / taumv(6)
              (hinf(6) - y(128)) / tauhv(6)
              KF * y(127) * y(128) - KB * y(129)
              KM * y(129) - KH * y(130)
              (actinf(6) - y(131)) / tauact(6)
              (inactinf(6) - y(132)) / tauinact(6)
              KFP * y(131) * y(132) - KBP * y(133)
              -Ilong(6) / capacL
            ];

% --------------------------------------------------------------------------

function plattau
     
global TAUM1 TAUM2 TAUMLIM STM TAUMDIF TAUH1 TAUHLIM TAUH2 STH FIGNUM AX1
global GAUSSTCL HY0 HA HOM HXC TAUNUM

for i = 1:71

    em(i) = i - 91;
    if (TAUMDIF ~= 0)
        taumv(i) = TAUM1 + TAUMDIF / (1 + exp((em(i) - TAUMLIM) * STM));
    else
        taumv(i) = TAUM1;
    end

    if (GAUSSTCL > 0)
        tauhv(i) = HY0 + HA  .* exp(-0.5 .* ((em(i)-HXC)./HOM).^2);
    else
        if (TAUH1 ~= TAUH2)
            tauhv(i) = TAUH1 + (TAUH2 - TAUH1) / (1 + exp((em(i) - TAUHLIM) * STH));
        else
            tauhv(i) = TAUH1;
        end
    end

end

% subplot(6,1,2);
figure (TAUNUM);
% hold on;
% ax2 = axes('Position',get(AX1,'Position'),...
%            'XAxisLocation','bottom',...
%            'YAxisLocation','right',...
%            'Color','none',...
%            'XColor','k','YColor','k');
% line(em,taumv,'color','k','linewidth',1,'parent',ax2);
% line(em,tauhv,'color','k','linewidth',1,'parent',ax2);
% axis([-90 -20 0 10]);
% ylabel('\tau_m & \tau_h (s)');
% set(gca,'YTick',[0 2 4 6 8 10],'ticklength',[.015,.015]);
plot(em,taumv,em,tauhv);
axis([-90 -20 0 10]);
figure(FIGNUM);
