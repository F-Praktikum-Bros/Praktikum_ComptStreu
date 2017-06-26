clc; close all; clear all; format long
g=@(E0) E0/0.511;  % E0 in MeV

% Daten einlesen und normieren
lt = importdata('zeiten'); % livetime

Ba = importdata('BA_ohneTarget.txt'); NBa=sum(Ba); errBa = sqrt(Ba);
Ba=Ba/lt(1); errBa=errBa/lt(1);

CsT = importdata('CS_mitTarget.txt'); NCsT=sum(CsT); errCsT = sqrt(CsT);
CsT = CsT/lt(2); errCsT=errCsT/lt(2);

Cs = importdata('CS_ohneTarget.txt'); NCs=sum(Cs); errCs = sqrt(Cs);
Cs=Cs/lt(3); errCs=errCs/lt(3);

Na = importdata('NA_ohneTarget.txt'); NNa=sum(Na); errNa = sqrt(Na);
Na=Na/lt(4); errNa=errNa/lt(4);

R = importdata('rauschen.txt'); NR=sum(R); errR = sqrt(R);
R=R/lt(5); errR=errR/lt(5);

dt=11; udt=1; % Todzeit in Mikrosekunden
% Realtimes
% time sind die hÃ¶heren zeiten, die zeiten in "zeiten" sind die niedrigeren - beachten welche reihenfolge in vektoren:
rt = importdata('Time'); rtBa=rt(1); rtCsT=rt(2); rtCs=rt(3); rtNa=rt(4); rtR=rt(5);
N=[NBa NCsT NCs NNa NR].';

diffTodzeit=(rt-lt)./N *10^6 -dt
size(Na)

% Rauschen korrigieren

Bac=Ba-R;   errBac = sqrt(errBa.^2 + errR.^2);
CsTc=CsT-R; errCsTc = sqrt(errCsT.^2 + errR.^2);
Csc=Cs-R;   errCsc = sqrt(errCs.^2 + errR.^2);
Nac=Na-R;   errNac = sqrt(errNa.^2 + errR.^2);

figure(26)
errorbar(1:length(Nac),Nac,errNac,'.','MarkerEdgeColor','black')


% Spektren unkorrigiert
figure(8)
errorbar(1:length(Ba),Ba,errBa,'.','MarkerEdgeColor','black')
xlabel('Kanal K'), ylabel('Anzahl Ereignisse pro Sekunde [1/s]')
xlim([0,length(R)])
saveas(gcf,'BaK.jpg')

figure(9)
errorbar(1:length(CsT),CsT,errCsT,'.','MarkerEdgeColor','black')
xlabel('Kanal K'), ylabel('Anzahl Ereignisse pro Sekunde [1/s]')
xlim([0,length(R)])
saveas(gcf,'CsTK.jpg')

figure(10)
errorbar(1:length(Cs),Cs,errCs,'.','MarkerEdgeColor','black')
xlabel('Kanal K'), ylabel('Anzahl Ereignisse pro Sekunde [1/s]')
xlim([0,length(R)])
saveas(gcf,'CsK.jpg')

figure(11)
errorbar(1:length(Nac),Nac,errNac,'.','MarkerEdgeColor','black')
xlabel('Kanal K'), ylabel('Anzahl Ereignisse pro Sekunde [1/s]')
xlim([0,length(R)])
saveas(gcf,'NaK.jpg')

% relevanter Bereich, Gauss-Peak
intan=210; intend=400;
Bar=Bac(250:330); errBar=errBac(235:315);
CsTr=CsTc(400:650); errCsTr=errCsTc(400:650);
Csr=Csc(400:650); errCsr=errCsc(400:650);
Nar=Nac(intan:intend); errNar=errNac(intan:intend);

sigmaumg=0.68;
% Plots der Gauss-Peaks
[fBar,gofBar] = fit(234+[1:length(Bar)].',Bar,'gauss1','Weights',errBar);
ci_fBar=confint(fBar,sigmaumg);
figure(1)
errorbar(234+[1:length(Bar)],Bar.',errBar,'.','MarkerEdgeColor','black')
hold on; plot(fBar);
ylabel('Anzahl Ereignisse pro Sekunde [1/s]'), xlabel('Kanal K'), xlim([233,318])
hold off
saveas(gcf,'peakBaK.jpg')

[fCsTr,gofCsTr] = fit(399+[1:length(CsTr)].',CsTr,'gauss1','Weights',errCsTr);
ci_fCsTr=confint(fCsTr,sigmaumg);
figure(2)
errorbar(399+[1:length(CsTr)],CsTr.',errCsTr,'.','MarkerEdgeColor','black')
hold on; plot(fCsTr);
ylabel('Anzahl Ereignisse pro Sekunde [1/s]'), xlabel('Kanal K'), xlim([387,647])
hold off

[fCsr,gofCsr] = fit(399+[1:length(Csr)].',Csr,'gauss1','Weights',errCsTr);
ci_fCsr=confint(fCsr,sigmaumg);
figure(3)
errorbar(399+[1:length(Csr)],Csr.',errCsr,'.','MarkerEdgeColor','black')
hold on; plot(fCsr);
ylabel('Anzahl Ereignisse pro Sekunde [1/s]'), xlabel('Kanal K'), xlim([387,647])
hold off
saveas(gcf,'peakCsrK.jpg')

[fNar,gofNar] = fit(intan-1+[1:length(Nar)].',Nar,'gauss1','Weights',errNar);
ci_fNar=confint(fNar,sigmaumg);
figure(4)
errorbar(intan-1+[1:length(Nar)],Nar.',errNar,'.','MarkerEdgeColor','black')
hold on; plot(fNar);
ylabel('Anzahl Ereignisse pro Sekunde [1/s]'), xlabel('Kanal K'), xlim([220,370])
hold off
saveas(gcf,'peakNarK.jpg')

% Peak-Positionen
bBar=fBar.b1; errbBar=0.5*(ci_fBar(2,2)-ci_fBar(1,2)); EBa=.356; % Kanal, MeV
bCsr=fCsr.b1; errbCsr=0.5*(ci_fCsr(2,2)-ci_fCsr(1,2)); ECs=.662; % Kanal, MeV
bNar=397.52; errbNar=0.278; ENa=.511; % Kanal, MeV

% Fit und Plot Kanal als Funktion von der Energie, K(E)
[flin,~]=fit([ENa EBa ECs].',[bNar bBar bCsr].','poly1','Weights',[errbNar errbBar errbCsr])
figure(5)
errorbar([EBa ENa  ECs],[bBar bNar  bCsr],[errbBar errbCsr errbNar],'.','MarkerEdgeColor','black')
xlabel('Energie E [MeV]'), ylabel('Kanal K')
hold on
plot(flin)
xlabel('Energie E [MeV]'), ylabel('Kanal K')
hold off
saveas(gcf,'eichkurveKvonE.jpg')
ci_flin=confint(flin,0.68);
up1=0.5*(ci_flin(2,1)-ci_flin(1,1));
up2=0.5*(ci_flin(2,2)-ci_flin(1,2));


% Energie als Funktion des Kanals, E(K), e1=1/p1, e2=-p2/p1
e1=1/flin.p1; ue1=0.5*(ci_flin(2,1)-ci_flin(1,1))/(flin.p1)^2;
e2=-flin.p2*e1; ue2=sqrt( (e2*e1^2 *0.5*(ci_flin(2,1)-ci_flin(1,1)))^2 + (0.5*(ci_flin(2,2)-ci_flin(1,2))*e1)^2 );
figure(6)
errorbar([EBa ENa  ECs],[bBar bNar  bCsr],[errbBar errbCsr errbNar],'.','MarkerEdgeColor','black')
set(gca,'CameraUpVector',[1,0,0],'YDir','reverse','XAxisLocation','top');
hold on
plot(e1*(bBar:0.01:bCsr)+e2,bBar:0.01:bCsr)
hold off
ylabel('Kanal K'), xlabel('Energie E [MeV]')
saveas(gcf,'eichkurveEvonK.jpg')
E=@(K) e1*K+e2; uE=@(K,uK) sqrt( (K*ue1).^2 + (ue2).^2 + (e1*uK).^2 );


% relative Aufloesung A(E)=h (E)/E des Szintillators, h(E)=2.35*sig(E)
sigBar  = fBar.c1;  usigBar   = 0.5*(ci_fBar(2,3)-ci_fBar(1,3));
sigCsTr = fCsTr.c1; usigCsTr  = 0.5*(ci_fCsTr(2,3)-ci_fCsTr(1,3));
sigCsr  = fCsr.c1;  usigCsr   = 0.5*(ci_fCsr(2,3)-ci_fCsr(1,3));
sigNar  = fNar.c1;  usigNar   = 0.5*(ci_fNar(2,3)-ci_fNar(1,3));

sigVec=[sigBar sigCsr sigNar]; usigVec= [usigBar usigCsr usigNar];
sigVec=sigVec;
usigVec=usigVec;
hVec=2.35*E(sigVec); uhVec=2.35*E(usigVec);
EVec=E([bBar bCsr bNar]);
A=hVec./EVec;
uEVec=uE([bBar bCsr bNar],[errbBar errbCsr errbNar]);
uA=sqrt( (uhVec./EVec).^2 + (uEVec.*hVec./(EVec).^2).^2 );

[xData, yData, weights] = prepareCurveData( EVec, A, uA );

% Set up fittype and options.
fitCsp = fittype( 'a/sqrt(x)+b', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.957506835434298 0.964888535199277];
opts.Weights = weights;

% Fit model to data.
[fitA, gof] = fit( xData, yData, fitCsp, opts );

figure(7)
errorbar(EVec,A,uA,'.','MarkerEdgeColor','black'); hold on
xlim([0.1,0.9]);
plot(fitA); hold off
xlabel('Energie E [MeV]'), ylabel('relative Auflösung A')
legend('Messpunkte \sigma(E)','Auflösungsfunktion A(E)=a\sqrt(E)+b','Interpreter','latex')
saveas(gcf,'relAufl.jpg')

disp('relative Auflösungen Ba Cs Na');
relA=@(E) (fitA.a)./sqrt(E)+fitA.b;

% Cs unkorrigiert Spektrum als Funktion der Energie (nicht mehr Kanal)
EVecCs=E(1:length(Cs));
figure(12)
errorbar(EVecCs,Cs,errCs,'.','MarkerEdgeColor','black')
xlim([0 EVecCs(end)+0.1]), xlabel('Energie E [MeV]'), ylabel('Counts pro Sekunde [1/s]')
saveas(gcf,'CsunkorrigiertE.jpg')

% Plot Cs mit und ohne Target (korrigiert)
figure(13)
errorbar(E(1:length(Csc)),Csc,errCsc,'.','MarkerEdgeColor','black')
hold on
%errorbar(E(1:length(CsTc)),CsTc,errCsTc,'x','MarkerEdgeColor','black')
xlim([0 1.5])
ylim([0 0.26])
set(gca,'XTick',[0:0.1:1.5]);
xlabel('Energie E [MeV]'), ylabel('Counts pro Sekunde [1/s]')
hold off
legend('Cs korrigiert','Location','northwest')
saveas(gcf,'Csk.jpg')


%Plot Cs ohne Target (unkorrigiert)
figure(26)
errorbar(E(1:length(Cs)),Cs,errCs,'.','MarkerEdgeColor','black')
xlim([0 1.5])
ylim([0 0.7])
set(gca,'XTick',[0:0.1:1.5]);
xlabel('Energie E [MeV]'), ylabel('Counts pro Sekunde [1/s]')
legend('Cs unkorrigiert','Location','northwest')
saveas(gcf,'Csuk.jpg')

% Plot Cs mit und ohne Target im Bereich des Peaks
Csp=Csc (305:680).'; errCsp=errCsc (305:680).'; CsTp=CsTc (305:680).'; errCsTp=errCsTc (305:680).';
EpVec=E(304)+E (1:length(Csp)).';
% figure(14)
% errorbar(EpVec,Csp,errCsp,'x','MarkerEdgeColor','black')
% hold on
% errorbar(EpVec,CsTp,errCsTp,'.','MarkerEdgeColor','black')
% hold off
% xlabel('Energie E [MeV]'), ylabel('Anzahl Ereignisse pro Sekunde [1/s]')

% Plot Rauschen
figure(17)
errorbar([1:length(R)],R,errR,'.','MarkerEdgeColor','black')
xlabel('Kanal K'), ylabel('Anzahl Ereignisse pro Sekunde [1/s]')
xlim([0,length(R)])
saveas(gcf,'Rauschen.jpg')


% Fit Compton WQ ohne Target
[fitCsp, ~] = createFitCsp(EpVec, Csp, errCsp)

% Fit Compton WQ mit Target
[fitCsTp, ~] = createFitCsTp(EpVec, CsTp, errCsTp)
figure(18)
hold on
errorbar(EpVec,Csp,errCsp,'.','MarkerEdgeColor','black')
errorbar(EpVec,CsTp,errCsTp,'x','MarkerEdgeColor','black')
plot(fitCsp)
plot(fitCsTp)
hold off
xlabel('Energie [MeV]'), ylabel('Anzahl Ereignisse pro Sekunde [1/s]')
xlim([E(328) E(697)])
legend('Cs + Fit','CsT + Fit','Location','northwest')
saveas(gcf,'FitsComptonCs.jpg')

ci_fitCsp=confint(fitCsp,0.68);
ci_fitCsTp=confint(fitCsTp,0.68);
uAoT=0.5*(ci_fitCsp(2,1)-ci_fitCsp(1,1));
uAmT=0.5*(ci_fitCsTp(2,1)-ci_fitCsTp(1,1));


% Compton-WQ
m= 4.4804*10^(-23); rho=2.72; Z= 13; x= 5; % cm ux=0.1cm
CWQ=m/(rho*Z*x)*log(fitCsp.a/fitCsTp.a)*10^24, % barn
uCWQ= sqrt( (m/(rho*Z*x^2)*log(fitCsp.a/fitCsTp.a)*0.1)^2 + (m/(rho*Z*x*fitCsp.a)*uAoT)^2 + (m/(rho*Z*x*fitCsTp.a)*uAmT)^2)*10^24

% Korrektur Raumwinkel
theta=atan(1/5); fun = @(t) (1./(1+g (ECs).*(1-cos(t)))).^2 .*( (1+g(ECs)*(1-cos(t))) + 1./(1+g (ECs).*(1-cos(t)))-sin (t).^2 ).*sin(t)
CWQKorr = pi*0.08*integral(fun,0,theta) % barn
CWQ = CWQ - CWQKorr;


% Plot totaler Wirkungsquerschnitt
sigma = @(g) 2*pi*8*10^(-30)* ((1+g)./g.^2 .* (2*(1+g)./(1+2*g) - 1./g .*log(1+2*g)) +1./(2*g) .*log(1+2*g) - (1+3*g)./(1+2*g).^2)*10^(28); % barn
E0Vec=0.1:0.01:100; % MeV
gVec=g(E0Vec);        % Vektor mit Gammawerten
sigmaVec=sigma(gVec);
figure(16)
plot(gVec,sigmaVec)
hold on
h = errorbar(g(ECs),CWQ,uCWQ,'.');
set(get(h,'Parent'),'XScale','log', 'YScale', 'log')
hold off
xlabel('Gamma'),ylabel('Wirkungsquerschnitt  [barn]')
saveas(gcf,'TotalerWQ.jpg')

% natuerliche Linienbreite: ~10^-24
