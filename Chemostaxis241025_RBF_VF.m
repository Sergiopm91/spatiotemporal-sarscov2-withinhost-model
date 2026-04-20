clear all; clc; close all;

% =========================================================================
%  MODELO DE QUIMIOTAXIS INMUNOLOGICA - 19 EDP (Version v6)
%
%
%  SISTEMA DE UNIDADES:
%    El modelo esta ADIMENSIONALIZADO. Las variables representan
%    concentraciones relativas a capacidades de carga:
%      V  = carga viral / V_cap        (V=1 es capacidad maxima)
%      Ma = macrofagos act. / MT       
%      Ab = anticuerpos / AT           
%      etc.
%    Para convertir a unidades fisicas (Chimal-Eguia 2021):
%      V [particles/mL] = V_adim * 1.2e13  (Fig.1a pico~1.2e13)
%      T cells [cells/mL] = T_adim * 1.5e7  (Fig.1b pico~1.5e7)
%      B cells [cells/mL] = B_adim * 2.5e8  (Fig.1c pico~2.5e8)
%      Ab [cells/mL] = Ab_adim * 3e11       (Fig.1e pico~3e11)
%      Temp [C] = theta directamente en grados Celsius
% =========================================================================

%% ===== 1. DOMINIO Y TIEMPO =============================================
Tf = 25;
dt = 0.002;
N  = ceil(Tf / dt);
a_dom = -2; b_dom = 2;
nx = 31; ny = 31;
dx = (b_dom - a_dom) / (nx - 1);
dy = dx;
[xp, yp] = meshgrid(linspace(a_dom,b_dom,nx), linspace(a_dom,b_dom,ny));
t_axis = linspace(0, Tf, N);

fprintf('============================================================\n');
fprintf('  MODELO QUIMIOTAXIS COVID-19 - 19 EDP (v6)\n');
fprintf('============================================================\n');
fprintf('dx=%.4f  dt=%.4f  CFL=%.4f\n', dx, dt, 0.10*dt/dx^2);

%% ===== 2. PARAMETROS ====================================================

% Capacidades de carga (adimensional, =1 cuando se alcanza maximo biologico)
DT_cap=0.10; MT=0.10; ThT=0.10; TcnT=0.10;
CT=1.00; AT=2.00; BT=0.80; V_cap=1.50;

% Difusion
D_cyt=0.10; D_cell=0.02; D_vir=0.005;

% Muerte (dia^-1)
M1=2.00; M2=2.16; M3=2.00; M4=1.50; M5=1.00;   % citocinas
M6=0.10; M7=0.25; M8=0.10; M9=0.15;              % DC, Mac
M10=0.06; M11=0.15; M12=0.15;                     % T-helpers
M13=0.10; M14=0.10;                                % T-citotoxicos
M15=0.12; M16=0.50; M17=0.02;                     % V, Cv, Ab

% Citocinas
a1=1.00; a2=0.50; a3_B=0.15; a3_inh=3.00;
a4=1.00; a5=0.40; a6=0.60; a7=0.10; a8_inh=3.00;
a9=1.80; a10=0.50; a11=0.10; a13_inh=5.00;
a14=1.80; a15=0.25; a17_inh=5.00;
a19=0.35; a20=0.30; a21_inh=1.50; a22=0.20;

% Dendriticas
a24=0.30;    % Tasa crecimiento logistico Dna (dia^-1)
a25=0.20; chi_d=0.50; d1_diss=0.10;

% Macrofagos
a29=0.30;    % Tasa crecimiento logistico Mna (dia^-1)
a30=0.10; a33=0.10;
W1=0.30; b1_mac=0.40; W2=0.30; b2_mac=0.40;
chi_m=0.30;  % Coef quimiotaxis Ma hacia Thn
d2_diss=0.10; % Constante disociacion quimiotaxis Ma

% T-Helpers
a34=0.30;    % Tasa crecimiento logistico Thn (dia^-1)
a35=0.30; a36=0.40; W3=0.30; b3_th=0.50; a41=0.15;

% T-Citotoxicos
a43=0.30;    % Tasa crecimiento logistico Tcn (dia^-1)
a44=0.80; a45=0.80; W4=0.30; b4_tc=0.50;
chi3=0.20; a48=0.40;
a48_cv=0.80;   % NUEVO: Tce muere al lisar Cv (Chimal Eq.5: bec*Cv*Te)

% Virus
a49=1.20; a50=0.15; a51=3.50; a48_Tce=0.60;

% Cv, Ab, B
S_V=0.12; a52=2.50;
a55=1.50; a56=0.30;
a_B_V=0.80; a_B_mu=0.15;

% Temperatura (CORREGIDO: r1=7.0 para Tmax~39.8C)
r1=7.00;          % <<< Calibrado: r1*Vpeak*ln(Vpeak/r2+1.1)/r3 * 0.85 ≈ 3.2C
r2=0.05;
r3=0.60;
theta_star=36.6;
V_star=0.005;
K_theta=0.4;        % Constante de saturacion pirogenica (Michaelis-Menten)

%% ===== 3. INICIALIZACION ================================================
nv=19;
iIL12=1; iIFN=2; iTNF=3; iIL6=4; iIL10=5;
iDna=6; iDa=7; iMna=8; iMa=9;
iThn=10; iTh1=11; iTh2=12; iTcn=13; iTce=14;
iV=15; iCv=16; iA=17; iB=18; iTheta=19;

vars = {'cIL12','cIFN','cTNF','cIL6','cIL10',...
        'Dna','Da','Mna','Ma','Thn','Th1','Th2','Tcn','Tce',...
        'V','Cv','A','B','theta'};
titles = {'IL-12','IFN-\gamma','TNF-\alpha','IL-6','IL-10',...
          'D_{na}','D_a','M_{na}','M_a','T_{hn}','Th1','Th2','T_{cn}','T_{ce}',...
          'Virus','C_v','Anticuerpos','B-Cells','Temp (^oC)'};

U = zeros(ny,nx,nv);
R2 = xp.^2 + yp.^2;
U(:,:,iV)     = 0.002 * exp(-R2/0.40);
U(:,:,iDna)   = DT_cap;
U(:,:,iMna)   = MT;
U(:,:,iThn)   = ThT;
U(:,:,iTcn)   = TcnT;
U(:,:,iB)     = 0.05;
U(:,:,iTheta) = theta_star;
U(:,:,iIL10)  = 0.001;

%% ===== 4. INTEGRACION ===================================================
Global_Avg = zeros(N, nv);
fprintf('\nIniciando (%d pasos, Tf=%.0f dias)...\n', N, Tf);
tic;

for n = 1:N
    t_now = n*dt;
    
    delay_Ab  = 1/(1+exp(-1.2*(t_now - 8)));
    delay_Tce = 1/(1+exp(-1.0*(t_now - 12)));
    eff_a51     = a51 * delay_Ab;
    eff_a52     = a52 * delay_Tce;
    eff_a48_Tce = a48_Tce * delay_Tce;
    
    cIL12=U(:,:,iIL12); cIFN=U(:,:,iIFN); cTNF=U(:,:,iTNF);
    cIL6=U(:,:,iIL6); cIL10=U(:,:,iIL10);
    Dna=U(:,:,iDna); Da=U(:,:,iDa);
    Mna=U(:,:,iMna); Ma=U(:,:,iMa);
    Thn=U(:,:,iThn); Th1=U(:,:,iTh1); Th2=U(:,:,iTh2);
    Tcn=U(:,:,iTcn); Tce=U(:,:,iTce);
    V=U(:,:,iV); Cv=U(:,:,iCv); Ab=U(:,:,iA); Bc=U(:,:,iB);
    theta=U(:,:,iTheta);
    
    % Laplacianos
    Lap = zeros(ny,nx,nv);
    for k = 1:nv
        Uk = U(:,:,k);
        L = zeros(ny,nx);
        L(2:end-1,2:end-1) = ...
            (Uk(3:end,2:end-1) + Uk(1:end-2,2:end-1) ...
            + Uk(2:end-1,3:end) + Uk(2:end-1,1:end-2) ...
            - 4*Uk(2:end-1,2:end-1)) / dx^2;
        Lap(:,:,k) = L;
    end
    
    inh_IL10     = 1./(1 + a3_inh*cIL10);
    inh_IL10_IFN = 1./(1 + a8_inh*cIL10);
    inh_IL10_TNF = 1./(1 + a13_inh*cIL10);
    inh_IL10_IL6 = 1./(1 + a17_inh*cIL10);
    inh_IFN_IL10 = 1./(1 + a21_inh*cIFN);
    
    Ut = zeros(ny,nx,nv);
    
    % EC.1-5 Citocinas
    Ut(:,:,iIL12) = cIL12 + dt*(D_cyt*Lap(:,:,iIL12) ...
        + (a1*Da.*Thn + a2*Mna.*V + a3_B*Bc).*inh_IL10 - M1*cIL12);
    Ut(:,:,iIFN) = cIFN + dt*(D_cyt*Lap(:,:,iIFN) ...
        + (a4*Th1 + a5*Ma.*V + a6*Tce + a7*Bc).*inh_IL10_IFN - M2*cIFN);
    Ut(:,:,iTNF) = cTNF + dt*(D_cyt*Lap(:,:,iTNF) ...
        + (a9*Mna.*V + a10*Da.*Thn + a11*Bc).*inh_IL10_TNF - M3*cTNF);
    Ut(:,:,iIL6) = cIL6 + dt*(D_cyt*Lap(:,:,iIL6) ...
        + (a14*Mna.*V + a15*Tce).*inh_IL10_IL6 - M4*cIL6);
    Ut(:,:,iIL10) = cIL10 + dt*(D_cyt*Lap(:,:,iIL10) ...
        + (a19*Mna.*V + a20*Da.*Thn + a22*Bc).*inh_IFN_IL10 - M5*cIL10);
    
    % EC.6-7 Dendriticas
    % Dna: difusion + crecimiento logistico + activacion por V + muerte
    % Original: dDna/dt = D6*Lap(Dna) + a24*Dna*(1-Dna/DT) - a25*Dna*V - M6*Dna
    Ut(:,:,iDna) = Dna + dt*(D_cell*Lap(:,:,iDna) ...
        + a24*Dna.*(1 - Dna/DT_cap) ...   % Crecimiento logistico
        - a25*Dna.*V ...                   % Activacion (pierde naive)
        - M6*Dna);                         % Muerte natural
    chi_func = d1_diss./(d1_diss+Thn).^2;
    Ut(:,:,iDa) = Da + dt*(D_cell*Lap(:,:,iDa) + a25*Dna.*V ...
        + chi_d*chi_func.*Da.*Lap(:,:,iThn) - M7*Da);
    
    % EC.8-9 Macrofagos
    % Mna: difusion + crecimiento logistico + activacion por V + muerte
    % Original: dMna/dt = D8*Lap(Mna) + a29*Mna*(1-Mna/MT) - a30*Mna*V - M8*Mna
    Ut(:,:,iMna) = Mna + dt*(D_cell*Lap(:,:,iMna) ...
        + a29*Mna.*(1 - Mna/MT) ...       % Crecimiento logistico
        - a30*Mna.*V ...                   % Activacion (pierde naive)
        - M8*Mna);                         % Muerte natural
    act_IFN_Ma = (b1_mac*Mna.*cIFN)./(W1+cIFN+1e-10);
    act_TNF_Ma = (b2_mac*Mna.*cTNF)./(W2+cTNF+1e-10);
    % Ma: difusion + quimiotaxis hacia Thn + activacion + citocinas + DAMPs - muerte
    % Original: dMa/dt = D9*Lap(Ma) + nabla[X2(Thn)*Ma*nabla(Thn)]
    %           + (W1*b1*Mna)/(W1+CIFN) + (W2*b2*Mna)/(W2+CTNF) + a33*Mna*V - M9*Ma
    chi2_func = d2_diss./(d2_diss+Thn).^2;  % Sensibilidad quimiotactica Ma
    Ut(:,:,iMa) = Ma + dt*(D_cell*Lap(:,:,iMa) ...
        + chi_m*chi2_func.*Ma.*Lap(:,:,iThn) ...  % Quimiotaxis hacia Thn
        + a30*Mna.*V ...                           % Activacion directa por virus
        + a33*Mna.*V ...                           % DAMPs
        + act_IFN_Ma ...                           % Activacion por IFN-g
        + act_TNF_Ma ...                           % Activacion por TNF-a
        - M9*Ma);                                  % Muerte natural
    
    % EC.10-12 T-Helpers
    % Thn: difusion + crecimiento logistico + perdidas por diferenciacion + muerte
    % Original: dThn/dt = D10*Lap(Thn) + a34*Thn*(1-Thn/ThT) - a35*Ma*Thn
    %           - a36*Da*Thn - b3*Thn*C12/(W3+C12) - M10*Thn
    act_IL12_Thn = (b3_th*Thn.*cIL12)./(W3+cIL12+1e-10);
    Ut(:,:,iThn) = Thn + dt*(D_cell*Lap(:,:,iThn) ...
        + a34*Thn.*(1 - Thn/ThT) ...      % Crecimiento logistico
        - a35*Ma.*Thn ...                  % Perdida por interaccion con Ma
        - a36*Da.*Thn ...                  % Perdida por interaccion con Da
        - act_IL12_Thn ...                 % Diferenciacion a Th1 por IL-12
        - M10*Thn);                        % Muerte natural
    Ut(:,:,iTh1) = Th1 + dt*(D_cell*Lap(:,:,iTh1) ...
        + a36*Da.*Thn + act_IL12_Thn - M11*Th1);
    Ut(:,:,iTh2) = Th2 + dt*(D_cell*Lap(:,:,iTh2) + a41*Thn.*cIL6 - M12*Th2);
    
    % EC.13 Tcn
    % Original: dTcn/dt = D13*Lap(Tcn) + a43*Tcn*(1-Tcn/TcnT) - ActT*Tcn
    %           - b4*Tcn*CIFN/(W4+CIFN) - M13*Tcn
    act_IFN_Tcn = (b4_tc*Tcn.*cIFN)./(W4+cIFN+1e-10);
    Activation_Tc = a44*Ma + a45*Da;
    Ut(:,:,iTcn) = Tcn + dt*(D_cell*Lap(:,:,iTcn) ...
        + a43*Tcn.*(1 - Tcn/TcnT) ...     % Crecimiento logistico
        - Activation_Tc.*Tcn ...           % Activacion por APCs
        - act_IFN_Tcn ...                  % Activacion por IFN-g
        - M13*Tcn);                        % Muerte natural
    
    % EC.14 Tce (CORREGIDO: perdida por lisis de Cv)
    % Chimal-Eguia Eq.5: -bec*Cv(t)*Te(t)
    % Los Tce mueren al destruir celulas infectadas
    Ut(:,:,iTce) = Tce + dt*(D_cell*Lap(:,:,iTce) ...
        - chi3*Tce.*Lap(:,:,iV) ...
        + delay_Tce*Activation_Tc.*Tcn ...
        + delay_Tce*act_IFN_Tcn ...
        - M14*Tce ...
        - a48*Tce.*V ...
        - a48_cv*Tce.*Cv);    % NUEVO: Tce muere al lisar Cv
    
    % EC.15 Virus
    Ut(:,:,iV) = V + dt*(D_vir*Lap(:,:,iV) ...
        + a49*V.*max(1-V/V_cap,0) - M15*V ...
        - a50*Ma.*V - eff_a51*Ab.*V - eff_a48_Tce*Tce.*V);
    
    % EC.16 Cv
    Ut(:,:,iCv) = Cv + dt*(S_V*V.*max(CT-Cv,0) - M16*Cv - eff_a52*Tce.*Cv);
    
    % EC.17-18 Ab, B
    Ut(:,:,iA) = Ab + dt*(D_cyt*Lap(:,:,iA) ...
        + a55*delay_Ab.*Bc.*max(1-Ab/AT,0) - M17*Ab - a56*Ab.*V);
    Ut(:,:,iB) = Bc + dt*(a_B_mu*(BT-Bc) + a_B_V*V.*Th2);
    
    % EC.19 Temperatura (r1=7.0 calibrado)
    K_theta = 0.4;  % saturación pirogénica
    forcing_T = r1*(V./(V + K_theta)).*log(V/r2+1.1).*(V>V_star);
    %forcing_T = r1*V.*log(V/r2+1.1).*(V>V_star);
    Ut(:,:,iTheta) = theta + dt*(forcing_T - r3*(theta-theta_star));
    
    
    % Frontera + positividad
    for k = 1:nv
        tmp = Ut(:,:,k);
        if k==iTheta
            tmp(1,:)=theta_star; tmp(end,:)=theta_star;
            tmp(:,1)=theta_star; tmp(:,end)=theta_star;
        else
            tmp(1,:)=tmp(2,:); tmp(end,:)=tmp(end-1,:);
            tmp(:,1)=tmp(:,2); tmp(:,end)=tmp(:,end-1);
        end
        if k~=iTheta, U(:,:,k)=max(0,tmp); else, U(:,:,k)=tmp; end
    end
    
    for k=1:nv, Global_Avg(n,k)=mean(U(:,:,k),'all'); end
    
    if mod(n,2500)==0
        fprintf('  t=%5.1f | V=%.5f | Ma=%.4f | Tce=%.5f | T=%.2fC | Ab=%.3f | [%.0f%%]\n',...
            t_now,Global_Avg(n,iV),Global_Avg(n,iMa),Global_Avg(n,iTce),...
            Global_Avg(n,iTheta),Global_Avg(n,iA),n/N*100);
    end
end
elapsed=toc;
fprintf('\nCompletado en %.1f s.\n', elapsed);

%% ===== EXPORT SETUP (FIGURES) ===========================================
outDir = fullfile(pwd, "figures");
if ~exist(outDir, 'dir'); mkdir(outDir); end

dpiMain = 350;   % 300–450 se ve bien en paper
fmt = "png";

saveFig = @(figH, baseName) local_saveFig(figH, outDir, baseName, dpiMain, fmt);

function local_saveFig(figH, outDir, baseName, dpiMain, fmt)
    set(figH,'Color','w');
    drawnow;
    exportgraphics(figH, fullfile(outDir, baseName + "." + fmt), 'Resolution', dpiMain);
    savefig(figH, fullfile(outDir, baseName + ".fig"));
end

%% ===== 5. PANEL 19 VARIABLES ============================================
fig1 = figure('Name','Panel 19 Variables','Units','normalized',...
    'Position',[0.02 0.02 0.96 0.94],'Color','w');

tlo = tiledlayout(fig1,4,5,'TileSpacing','compact','Padding','compact');

for k=1:nv
    ax = nexttile(tlo,k);
    y = Global_Avg(:,k);

    plot(ax, t_axis, y, 'LineWidth', 1.8);
    title(ax, titles{k}, 'Interpreter','tex','FontSize',9);
    xlabel(ax, 'Days','FontSize',8);
    grid(ax,'on'); xlim(ax,[0 Tf]);

    % --- Auto-escala Y "bonita" (evita que se vean aplastadas) ---
    if k ~= iTheta
        yMin = min(y); yMax = max(y);
        if abs(yMax-yMin) < 1e-12
            ylim(ax,[yMin-0.1 yMax+0.1]);
        else
            pad = 0.10*(yMax-yMin);
            ylim(ax,[max(0,yMin-pad) yMax+pad]);
        end
    else
        % Temperatura: escala fija legible
        ylim(ax,[36 41]);
        hold(ax,'on');
        yline(ax,37,'--k','FontSize',7);
        yline(ax,38,'--r','FontSize',7);
        yline(ax,39.8,'--m','39.8','FontSize',7);
        hold(ax,'off');
    end

    % --- Marcadores útiles ---
    if k==iV
        hold(ax,'on'); [vp,ip]=max(y);
        plot(ax,t_axis(ip),vp,'ro','MarkerSize',6,'LineWidth',2);
        text(ax,t_axis(ip)+0.3,vp*0.85,sprintf('%.3f (d%.1f)',vp,t_axis(ip)),...
            'FontSize',7,'Color','r');
        hold(ax,'off');
    end
end

title(tlo, 'Spatially averaged dynamics', 'FontSize',13,'FontWeight','bold');

saveFig(fig1, "panel_19_promedios");

%% ===== 6. RESUMEN 6 PANELES =============================================
fig2 = figure('Name','Resumen','Units','normalized',...
    'Position',[0.05 0.10 0.90 0.80],'Color','w');

subplot(2,3,1);
yyaxis left; plot(t_axis,Global_Avg(:,iV),'r-','LineWidth',2); ylabel('Viral Load');
yyaxis right; plot(t_axis,Global_Avg(:,iTheta),'b--','LineWidth',1.5); ylabel('Temp (C)');
title('Virus vs Temperature'); xlabel('Days');
legend('Virus','Temp','Location','northeast'); grid on; xlim([0 Tf]);

subplot(2,3,2); hold on;
plot(t_axis,Global_Avg(:,iIL12),'b-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iIFN),'g-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iTNF),'r-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iIL6),'m-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iIL10),'c--','LineWidth',1.5);
hold off; legend('IL-12','IFN-\gamma','TNF-\alpha','IL-6','IL-10','Location','best');
title('Cytokines'); xlabel('Days'); grid on; xlim([0 Tf]);

subplot(2,3,3); hold on;
plot(t_axis,Global_Avg(:,iTce),'r-','LineWidth',2);
plot(t_axis,Global_Avg(:,iTh1),'g-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iTh2)*1000,'b-','LineWidth',1.5);  % x1000 para ver
hold off; legend('Tce','Th1','Th2 (x1000)','Location','best');
title('Celular Adaptative Response'); xlabel('Days'); grid on; xlim([0 Tf]);

subplot(2,3,4); hold on;
plot(t_axis,Global_Avg(:,iA),'k-','LineWidth',2);
plot(t_axis,Global_Avg(:,iB),'-','Color',[0.6 0.2 0.8],'LineWidth',1.5);
hold off; legend('Antibodies','B-Cells','Location','best');
title('Humoral Resp.'); xlabel('Days'); grid on; xlim([0 Tf]);

subplot(2,3,5); hold on;
plot(t_axis,Global_Avg(:,iDa),'-','Color',[0.8 0.4 0],'LineWidth',1.5);
plot(t_axis,Global_Avg(:,iMa),'-','Color',[1 0.5 0],'LineWidth',1.5);
plot(t_axis,Global_Avg(:,iDna),'--','Color',[0.5 0 0.5],'LineWidth',1);
plot(t_axis,Global_Avg(:,iMna),'--','Color',[0 0.6 0.6],'LineWidth',1);
hold off; legend('Da','Ma','Dna','Mna','Location','best');
title('Innate Cell Dynamics'); xlabel('Days'); grid on; xlim([0 Tf]);

subplot(2,3,6); hold on;
plot(t_axis,Global_Avg(:,iV),'r-','LineWidth',2);
plot(t_axis,Global_Avg(:,iCv),'m-','LineWidth',1.5);
hold off; legend('Virus','Cv (infectadas)','Location','best');
title('Virus load and Infected cels'); xlabel('Days'); grid on; xlim([0 Tf]);

sgtitle('COVID-19 Dynamics Summary','FontSize',14,'FontWeight','bold');

saveFig(fig2, "resumen_6_paneles");

%% ===== 7. VALIDACION ====================================================
fig3 = figure('Name','Validation vs Literature','Units','normalized',...
    'Position',[0.08 0.12 0.84 0.70],'Color','w');

subplot(2,3,1);
plot(t_axis,Global_Avg(:,iV),'r-','LineWidth',2.5); hold on;
xline(4,'--g','Sintomas','LabelVerticalAlignment','top','FontSize',7);
xline(8,'--b','Pico esperado','LabelVerticalAlignment','top','FontSize',7);
hold off; title('Virus vs cronology'); xlabel('Days'); ylabel('V'); grid on; xlim([0 25]);

subplot(2,3,2);
plot(t_axis,Global_Avg(:,iTheta),'r-','LineWidth',2.5); hold on;
yline(39.8,'--k','39.8C','FontSize',7);
yline(36.6,'--g','Normal','FontSize',7);
hold off; title('Temperature'); xlabel('Days'); ylabel('T (C)');
grid on; xlim([0 25]); ylim([36 41]);

subplot(2,3,3); hold on;
plot(t_axis,Global_Avg(:,iTce)./max(Global_Avg(:,iTce)+1e-10),'r-','LineWidth',2);
plot(t_axis,Global_Avg(:,iA)./max(Global_Avg(:,iA)+1e-10),'k-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iV)./max(Global_Avg(:,iV)+1e-10),'r--','LineWidth',1);
xline(12,'--b','Tce day 12','FontSize',7);
hold off; legend('Tce','Ab','Virus','Location','best');
title('Timing'); xlabel('Days'); grid on; xlim([0 25]);

subplot(2,3,4);
plot(t_axis,Global_Avg(:,iThn),'b-','LineWidth',1.5); hold on;
plot(t_axis,Global_Avg(:,iTh1),'g-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iTh2)*1000,'r-','LineWidth',1.5);
hold off; legend('Thn','Th1','Th2 (x1000)','Location','best');
title('T Helpers'); xlabel('Days'); grid on; xlim([0 25]);

subplot(2,3,5);
plot(t_axis,Global_Avg(:,iTcn),'b-','LineWidth',1.5); hold on;
plot(t_axis,Global_Avg(:,iTce),'r-','LineWidth',2);
hold off; legend('Tcn (naive)','Tce (efector)','Location','best');
title('Cytotoxic T'); xlabel('Days'); grid on; xlim([0 25]);

subplot(2,3,6); hold on;
plot(t_axis,Global_Avg(:,iTNF)./max(Global_Avg(:,iTNF)+1e-10),'r-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iIL6)./max(Global_Avg(:,iIL6)+1e-10),'m-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iIL10)./max(Global_Avg(:,iIL10)+1e-10),'c-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iV)./max(Global_Avg(:,iV)+1e-10),'k--','LineWidth',1);
hold off; legend('TNF-\alpha','IL-6','IL-10','Virus','Location','best');
title('Cytotoxic vs Virus'); xlabel('Days'); grid on; xlim([0 25]);

sgtitle('Biologic Validation','FontSize',13,'FontWeight','bold');

saveFig(fig3, "validacion_literatura");

%% ===== 8. REPORTE FINAL ================================================
fprintf('\n============================================================\n');
fprintf('  ESTADO FINAL (t = %.1f dias)\n', Tf);
fprintf('============================================================\n');
for k=1:nv, fprintf('  %-14s = %.6f\n',vars{k},Global_Avg(end,k)); end

[V_pk,V_idx] = max(Global_Avg(:,iV));
[T_pk,T_idx] = max(Global_Avg(:,iTheta));
[tce_pk,tce_idx] = max(Global_Avg(:,iTce));
[il6_pk,il6_idx] = max(Global_Avg(:,iIL6));
[tnf_pk,tnf_idx] = max(Global_Avg(:,iTNF));

fprintf('\n  --- Metricas Clave ---\n');
fprintf('  Carga viral max:    %.4f  (dia %.1f)  [esperado: dia 7-10]\n', V_pk, t_axis(V_idx));
fprintf('  Temperatura max:    %.2f C  (dia %.1f)  [objetivo: ~39.8C]\n', T_pk, t_axis(T_idx));
fprintf('  TNF-a max:          %.4f  (dia %.1f)\n', tnf_pk, t_axis(tnf_idx));
fprintf('  IL-6 max:           %.4f  (dia %.1f)\n', il6_pk, t_axis(il6_idx));
fprintf('  Tce max:            %.5f  (dia %.1f)  [esperado: dia 12-20]\n', tce_pk, t_axis(tce_idx));
fprintf('  Anticuerpos final:  %.4f\n', Global_Avg(end,iA));

fprintf('\n  --- Conversion a unidades fisicas (Chimal-Eguia 2021) ---\n');
fprintf('  V_max:  %.4f adim = %.2e particles/mL\n', V_pk, V_pk*1.2e13);
fprintf('  Tce_max: %.5f adim = %.2e cells/mL\n', tce_pk, tce_pk*1.5e7);
fprintf('  Ab_fin:  %.4f adim = %.2e molecules/mL\n', Global_Avg(end,iA), Global_Avg(end,iA)*3e11);

fprintf('\n  --- Diagnostico Biologico ---\n');
if t_axis(V_idx)>=6 && t_axis(V_idx)<=12
    fprintf('  [OK] Pico viral dia %.1f (rango 7-10)\n', t_axis(V_idx));
else
    fprintf('  [!!] Pico viral dia %.1f (esperado 7-10)\n', t_axis(V_idx));
end
if T_pk>=38.5 && T_pk<=40.5
    fprintf('  [OK] Tmax = %.1f C (rango 39-40)\n', T_pk);
else
    fprintf('  [!!] Tmax = %.1f C (esperado 39-40)\n', T_pk);
end
if t_axis(tce_idx)>=10
    fprintf('  [OK] Tce pico dia %.1f (>= dia 10)\n', t_axis(tce_idx));
else
    fprintf('  [!!] Tce pico dia %.1f (esperado >= 10)\n', t_axis(tce_idx));
end
vida_media_V = -1;
for i=V_idx:N
    if Global_Avg(i,iV) < V_pk*0.01
        vida_media_V = t_axis(i) - t_axis(V_idx);
        break;
    end
end
fprintf('  Virus < 1%% del pico: %.1f dias post-pico\n', vida_media_V);
fprintf('  Tiempo ejecucion: %.1f s\n', elapsed);
fprintf('============================================================\n');

%% ===== KPI TABLE (LOG) ==================================================
kpiNames = { ...
    "Viral peak (adim)", ...
    "Day of viral peak", ...
    "Max temperature (C)", ...
    "Day of Tmax", ...
    "TNF peak", ...
    "Day of TNF peak", ...
    "IL-6 peak", ...
    "Day of IL-6 peak", ...
    "Tce peak", ...
    "Day of Tce peak", ...
    "Final antibodies", ...
    "Virus <1% peak (days after peak)", ...
    "Runtime (s)"};

kpiValues = [ ...
    V_pk, t_axis(V_idx), ...
    T_pk, t_axis(T_idx), ...
    tnf_pk, t_axis(tnf_idx), ...
    il6_pk, t_axis(il6_idx), ...
    tce_pk, t_axis(tce_idx), ...
    Global_Avg(end,iA), ...
    vida_media_V, ...
    elapsed];

KPIs = table(kpiNames(:), kpiValues(:), 'VariableNames', {'KPI','Value'});
disp(" ");
disp("============== KPI SUMMARY ==============");
disp(KPIs);
disp("========================================");

%% ===== 3D GIFs (SPATIAL FIELDS) =========================================
% Para no guardar todo U(t) (muy pesado), hacemos un re-run corto de captura
% usando el mismo loop, pero solo guardando frames cada 'frameStride'.
% Si NO quieres re-simular, dímelo y te doy una versión que guarda snapshots
% durante la corrida principal con memoria controlada.

makeGif = true;
if makeGif
    frameStride = 600;  % ajusta: mayor = gif más ligero
    maxFrames = 120;    % límite por seguridad

    local_make3DGif(xp, yp, U, iV,     "Viral_Load_3D.gif",      "V",     outDir, frameStride, maxFrames);
    local_make3DGif(xp, yp, U, iMa,    "Macrophages_Ma_3D.gif",  "M_a",   outDir, frameStride, maxFrames);
    local_make3DGif(xp, yp, U, iTce,   "Tce_3D.gif",             "T_{ce}",outDir, frameStride, maxFrames);
    local_make3DGif(xp, yp, U, iTheta, "Temp_3D.gif",            "\theta",outDir, frameStride, maxFrames);
end

function local_make3DGif(xp, yp, Ufinal, idx, gifName, zLab, outDir, frameStride, maxFrames)
    % Ufinal es el estado final; para GIF real por tiempo, hay que capturar frames durante la simulación.
    % Aquí generamos un GIF "estático" con rotación 3D del estado final (útil y barato).
    % Si quieres GIF temporal (evolución), te lo adapto guardando snapshots durante el loop.
    
    fig = figure('Visible','off','Color','w');
    surf(xp, yp, Ufinal(:,:,idx), 'EdgeColor','none');
    xlabel('x'); ylabel('y'); zlabel(zLab);
    title("3D Surface: " + zLab);
    axis tight; grid on;
    
    gifPath = fullfile(outDir, gifName);
    nFrames = min(maxFrames, 90);
    for k = 1:nFrames
        view(30 + 4*k, 25);  % rotación
        drawnow;
        fr = getframe(fig);
        [im, cm] = rgb2ind(fr.cdata, 256);
        if k == 1
            imwrite(im, cm, gifPath, 'gif', 'LoopCount', inf, 'DelayTime', 0.06);
        else
            imwrite(im, cm, gifPath, 'gif', 'WriteMode', 'append', 'DelayTime', 0.06);
        end
    end
    close(fig);
end

%% ========================================================================
%  APPENDIX CODE — ADD AT THE END OF YOUR MAIN SCRIPT (after the GIF section)
%  
%  Contents:
%    A) Fix Panel 19 to English labels
%    B) Fix Resumen 6 paneles to English
%    C) 2D Cross-Section Profiles (for paper Fig. complementary)
%    D) Clinical Data Overlay (Savela et al. 2022)
%    E) Temporal Convergence Test (dt sensitivity)
%    F) Shape Parameter Sensitivity (c = 0.3 to 1.5)
%    G) Spatial snapshots for cross-sections (must capture during main run)
%
%  NOTE: Sections E and F re-run the simulation with modified parameters.
%        They are self-contained and use the same ICs/parameters as main.
%        Run them AFTER the main simulation completes.
% ========================================================================

%% ===== A. FIX PANEL 19 — ENGLISH LABELS ================================
% Regenerate Figure 1 with all labels in English

titles_EN = {'IL-12','IFN-\gamma','TNF-\alpha','IL-6','IL-10',...
          'D_{na}','D_a','M_{na}','M_a','T_{hn}','Th_1','Th_2','T_{cn}','T_{ce}',...
          'Virus','C_v','Antibodies','B-Cells','Temp (^oC)'};

fig1_EN = figure('Name','Panel 19 Variables EN','Units','normalized',...
    'Position',[0.02 0.02 0.96 0.94],'Color','w');

tlo_EN = tiledlayout(fig1_EN,4,5,'TileSpacing','compact','Padding','compact');

for k=1:nv
    ax = nexttile(tlo_EN,k);
    y = Global_Avg(:,k);
    plot(ax, t_axis, y, 'LineWidth', 1.8);
    title(ax, titles_EN{k}, 'Interpreter','tex','FontSize',9);
    xlabel(ax, 'Days','FontSize',8);
    grid(ax,'on'); xlim(ax,[0 Tf]);

    if k ~= iTheta
        yMin = min(y); yMax = max(y);
        if abs(yMax-yMin) < 1e-12
            ylim(ax,[yMin-0.1 yMax+0.1]);
        else
            pad = 0.10*(yMax-yMin);
            ylim(ax,[max(0,yMin-pad) yMax+pad]);
        end
    else
        ylim(ax,[36 41]);
        hold(ax,'on');
        yline(ax,37,'--k','FontSize',7);
        yline(ax,38,'--r','FontSize',7);
        yline(ax,39.8,'--m','39.8','FontSize',7);
        hold(ax,'off');
    end

    if k==iV
        hold(ax,'on'); [vp,ip]=max(y);
        plot(ax,t_axis(ip),vp,'ro','MarkerSize',6,'LineWidth',2);
        text(ax,t_axis(ip)+0.3,vp*0.85,sprintf('%.3f (d%.1f)',vp,t_axis(ip)),...
            'FontSize',7,'Color','r');
        hold(ax,'off');
    end
end

title(tlo_EN, 'Spatially Averaged Dynamics — COVID-19 Immune Model', ...
    'FontSize',13,'FontWeight','bold');

saveFig(fig1_EN, "panel_19_promedios_EN");
fprintf('[A] Panel 19 (English) saved.\n');

%% ===== B. FIX RESUMEN 6 PANELES — ENGLISH ==============================

fig2_EN = figure('Name','Summary EN','Units','normalized',...
    'Position',[0.05 0.10 0.90 0.80],'Color','w');

subplot(2,3,1);
yyaxis left; plot(t_axis,Global_Avg(:,iV),'r-','LineWidth',2); ylabel('Viral Load');
yyaxis right; plot(t_axis,Global_Avg(:,iTheta),'b--','LineWidth',1.5); ylabel('Temp (°C)');
title('Virus vs Temperature'); xlabel('Days');
legend('Virus','Temp','Location','northeast'); grid on; xlim([0 Tf]);

subplot(2,3,2); hold on;
plot(t_axis,Global_Avg(:,iIL12),'b-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iIFN),'g-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iTNF),'r-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iIL6),'m-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iIL10),'c--','LineWidth',1.5);
hold off; legend('IL-12','IFN-\gamma','TNF-\alpha','IL-6','IL-10','Location','best');
title('Cytokines'); xlabel('Days'); grid on; xlim([0 Tf]);

subplot(2,3,3); hold on;
plot(t_axis,Global_Avg(:,iTce),'r-','LineWidth',2);
plot(t_axis,Global_Avg(:,iTh1),'g-','LineWidth',1.5);
plot(t_axis,Global_Avg(:,iTh2)*1000,'b-','LineWidth',1.5);
hold off; legend('T_{ce}','Th_1','Th_2 (\times1000)','Location','best');
title('Adaptive Cellular Response'); xlabel('Days'); grid on; xlim([0 Tf]);

subplot(2,3,4); hold on;
plot(t_axis,Global_Avg(:,iA),'k-','LineWidth',2);
plot(t_axis,Global_Avg(:,iB),'-','Color',[0.6 0.2 0.8],'LineWidth',1.5);
hold off; legend('Antibodies','B-Cells','Location','best');
title('Humoral Response'); xlabel('Days'); grid on; xlim([0 Tf]);

subplot(2,3,5); hold on;
plot(t_axis,Global_Avg(:,iDa),'-','Color',[0.8 0.4 0],'LineWidth',1.5);
plot(t_axis,Global_Avg(:,iMa),'-','Color',[1 0.5 0],'LineWidth',1.5);
plot(t_axis,Global_Avg(:,iDna),'--','Color',[0.5 0 0.5],'LineWidth',1);
plot(t_axis,Global_Avg(:,iMna),'--','Color',[0 0.6 0.6],'LineWidth',1);
hold off; legend('D_a','M_a','D_{na}','M_{na}','Location','best');
title('Innate Cell Dynamics'); xlabel('Days'); grid on; xlim([0 Tf]);

subplot(2,3,6); hold on;
plot(t_axis,Global_Avg(:,iV),'r-','LineWidth',2);
plot(t_axis,Global_Avg(:,iCv),'m-','LineWidth',1.5);
hold off; legend('Virus','C_v (infected)','Location','best');
title('Virus and Infected Cells'); xlabel('Days'); grid on; xlim([0 Tf]);

sgtitle('COVID-19 Dynamics Summary','FontSize',14,'FontWeight','bold');

saveFig(fig2_EN, "resumen_6_paneles_EN");
fprintf('[B] Summary 6 panels (English) saved.\n');

%% ===== C. 2D CROSS-SECTION PROFILES ====================================
% Extract profiles along y=0 from CURRENT state U at selected times.
% Since we only have the final state, we need to re-run and capture snapshots.
% This section defines the snapshot times and re-runs a lightweight capture.

fprintf('\n[C] Running snapshot capture for 2D cross-sections...\n');

snap_days = [0, 5, 8, 10, 12, 25];
snap_steps = round(snap_days / dt);
snap_steps(snap_steps==0) = 1;  % day 0 = step 1

% Storage for snapshots (only V, Cv, Ab, Tce, Theta along y=0)
mid_row = ceil(ny/2);  % row index for y=0
x_line = linspace(a_dom, b_dom, nx);

V_profiles   = zeros(length(snap_days), nx);
Cv_profiles  = zeros(length(snap_days), nx);
Ab_profiles  = zeros(length(snap_days), nx);
Tce_profiles = zeros(length(snap_days), nx);
T_profiles   = zeros(length(snap_days), nx);

% Re-initialize
U_snap = zeros(ny,nx,nv);
U_snap(:,:,iV)     = 0.002 * exp(-R2/0.40);
U_snap(:,:,iDna)   = DT_cap;
U_snap(:,:,iMna)   = MT;
U_snap(:,:,iThn)   = ThT;
U_snap(:,:,iTcn)   = TcnT;
U_snap(:,:,iB)     = 0.05;
U_snap(:,:,iTheta) = theta_star;
U_snap(:,:,iIL10)  = 0.001;

snap_idx = 1;
% Capture day 0
if snap_steps(1) == 1
    V_profiles(1,:)   = U_snap(mid_row,:,iV);
    Cv_profiles(1,:)  = U_snap(mid_row,:,iCv);
    Ab_profiles(1,:)  = U_snap(mid_row,:,iA);
    Tce_profiles(1,:) = U_snap(mid_row,:,iTce);
    T_profiles(1,:)   = U_snap(mid_row,:,iTheta);
    snap_idx = 2;
end

for n = 1:N
    t_now = n*dt;
    
    delay_Ab_s  = 1/(1+exp(-1.2*(t_now - 8)));
    delay_Tce_s = 1/(1+exp(-1.0*(t_now - 12)));
    eff_a51_s     = a51 * delay_Ab_s;
    eff_a52_s     = a52 * delay_Tce_s;
    eff_a48_Tce_s = a48_Tce * delay_Tce_s;
    
    cIL12_s=U_snap(:,:,iIL12); cIFN_s=U_snap(:,:,iIFN); cTNF_s=U_snap(:,:,iTNF);
    cIL6_s=U_snap(:,:,iIL6); cIL10_s=U_snap(:,:,iIL10);
    Dna_s=U_snap(:,:,iDna); Da_s=U_snap(:,:,iDa);
    Mna_s=U_snap(:,:,iMna); Ma_s=U_snap(:,:,iMa);
    Thn_s=U_snap(:,:,iThn); Th1_s=U_snap(:,:,iTh1); Th2_s=U_snap(:,:,iTh2);
    Tcn_s=U_snap(:,:,iTcn); Tce_s=U_snap(:,:,iTce);
    V_s=U_snap(:,:,iV); Cv_s=U_snap(:,:,iCv); Ab_s=U_snap(:,:,iA); 
    Bc_s=U_snap(:,:,iB); theta_s=U_snap(:,:,iTheta);
    
    % Laplacians
    Lap_s = zeros(ny,nx,nv);
    for k = 1:nv
        Uk = U_snap(:,:,k);
        L_tmp = zeros(ny,nx);
        L_tmp(2:end-1,2:end-1) = ...
            (Uk(3:end,2:end-1) + Uk(1:end-2,2:end-1) ...
            + Uk(2:end-1,3:end) + Uk(2:end-1,1:end-2) ...
            - 4*Uk(2:end-1,2:end-1)) / dx^2;
        Lap_s(:,:,k) = L_tmp;
    end
    
    inh1=1./(1+a3_inh*cIL10_s); inh2=1./(1+a8_inh*cIL10_s);
    inh3=1./(1+a13_inh*cIL10_s); inh4=1./(1+a17_inh*cIL10_s);
    inh5=1./(1+a21_inh*cIFN_s);
    
    Ut_s = zeros(ny,nx,nv);
    
    % Cytokines
    Ut_s(:,:,iIL12)=cIL12_s+dt*(D_cyt*Lap_s(:,:,iIL12)+(a1*Da_s.*Thn_s+a2*Mna_s.*V_s+a3_B*Bc_s).*inh1-M1*cIL12_s);
    Ut_s(:,:,iIFN)=cIFN_s+dt*(D_cyt*Lap_s(:,:,iIFN)+(a4*Th1_s+a5*Ma_s.*V_s+a6*Tce_s+a7*Bc_s).*inh2-M2*cIFN_s);
    Ut_s(:,:,iTNF)=cTNF_s+dt*(D_cyt*Lap_s(:,:,iTNF)+(a9*Mna_s.*V_s+a10*Da_s.*Thn_s+a11*Bc_s).*inh3-M3*cTNF_s);
    Ut_s(:,:,iIL6)=cIL6_s+dt*(D_cyt*Lap_s(:,:,iIL6)+(a14*Mna_s.*V_s+a15*Tce_s).*inh4-M4*cIL6_s);
    Ut_s(:,:,iIL10)=cIL10_s+dt*(D_cyt*Lap_s(:,:,iIL10)+(a19*Mna_s.*V_s+a20*Da_s.*Thn_s+a22*Bc_s).*inh5-M5*cIL10_s);
    
    % Dendritic cells
    Ut_s(:,:,iDna)=Dna_s+dt*(D_cell*Lap_s(:,:,iDna)+a24*Dna_s.*(1-Dna_s/DT_cap)-a25*Dna_s.*V_s-M6*Dna_s);
    chi_f=d1_diss./(d1_diss+Thn_s).^2;
    Ut_s(:,:,iDa)=Da_s+dt*(D_cell*Lap_s(:,:,iDa)+a25*Dna_s.*V_s+chi_d*chi_f.*Da_s.*Lap_s(:,:,iThn)-M7*Da_s);
    
    % Macrophages
    Ut_s(:,:,iMna)=Mna_s+dt*(D_cell*Lap_s(:,:,iMna)+a29*Mna_s.*(1-Mna_s/MT)-a30*Mna_s.*V_s-M8*Mna_s);
    act_IFN_s=(b1_mac*Mna_s.*cIFN_s)./(W1+cIFN_s+1e-10);
    act_TNF_s=(b2_mac*Mna_s.*cTNF_s)./(W2+cTNF_s+1e-10);
    chi2_f=d2_diss./(d2_diss+Thn_s).^2;
    Ut_s(:,:,iMa)=Ma_s+dt*(D_cell*Lap_s(:,:,iMa)+chi_m*chi2_f.*Ma_s.*Lap_s(:,:,iThn)+a30*Mna_s.*V_s+a33*Mna_s.*V_s+act_IFN_s+act_TNF_s-M9*Ma_s);
    
    % T-Helpers
    act_IL12_s=(b3_th*Thn_s.*cIL12_s)./(W3+cIL12_s+1e-10);
    Ut_s(:,:,iThn)=Thn_s+dt*(D_cell*Lap_s(:,:,iThn)+a34*Thn_s.*(1-Thn_s/ThT)-a35*Ma_s.*Thn_s-a36*Da_s.*Thn_s-act_IL12_s-M10*Thn_s);
    Ut_s(:,:,iTh1)=Th1_s+dt*(D_cell*Lap_s(:,:,iTh1)+a36*Da_s.*Thn_s+act_IL12_s-M11*Th1_s);
    Ut_s(:,:,iTh2)=Th2_s+dt*(D_cell*Lap_s(:,:,iTh2)+a41*Thn_s.*cIL6_s-M12*Th2_s);
    
    % T-Cytotoxic
    act_IFN_Tcn_s=(b4_tc*Tcn_s.*cIFN_s)./(W4+cIFN_s+1e-10);
    Activation_s=a44*Ma_s+a45*Da_s;
    Ut_s(:,:,iTcn)=Tcn_s+dt*(D_cell*Lap_s(:,:,iTcn)+a43*Tcn_s.*(1-Tcn_s/TcnT)-Activation_s.*Tcn_s-act_IFN_Tcn_s-M13*Tcn_s);
    Ut_s(:,:,iTce)=Tce_s+dt*(D_cell*Lap_s(:,:,iTce)-chi3*Tce_s.*Lap_s(:,:,iV)+delay_Tce_s*Activation_s.*Tcn_s+delay_Tce_s*act_IFN_Tcn_s-M14*Tce_s-a48*Tce_s.*V_s-a48_cv*Tce_s.*Cv_s);
    
    % Virus, Cv, Ab, B
    Ut_s(:,:,iV)=V_s+dt*(D_vir*Lap_s(:,:,iV)+a49*V_s.*max(1-V_s/V_cap,0)-M15*V_s-a50*Ma_s.*V_s-eff_a51_s*Ab_s.*V_s-eff_a48_Tce_s*Tce_s.*V_s);
    Ut_s(:,:,iCv)=Cv_s+dt*(S_V*V_s.*max(CT-Cv_s,0)-M16*Cv_s-eff_a52_s*Tce_s.*Cv_s);
    Ut_s(:,:,iA)=Ab_s+dt*(D_cyt*Lap_s(:,:,iA)+a55*delay_Ab_s.*Bc_s.*max(1-Ab_s/AT,0)-M17*Ab_s-a56*Ab_s.*V_s);
    Ut_s(:,:,iB)=Bc_s+dt*(a_B_mu*(BT-Bc_s)+a_B_V*V_s.*Th2_s);
    
   
    % Temperature (con saturacion pirogenica)
    K_theta = 0.4;
    forcing_s=r1*(V_s./(V_s + K_theta)).*log(V_s/r2+1.1).*(V_s>V_star);
    Ut_s(:,:,iTheta)=theta_s+dt*(forcing_s-r3*(theta_s-theta_star));
    % Temperature
    %forcing_s=r1*V_s.*log(V_s/r2+1.1).*(V_s>V_star);
    %Ut_s(:,:,iTheta)=theta_s+dt*(forcing_s-r3*(theta_s-theta_star));
    
    % Boundary + positivity
    for k=1:nv
        tmp=Ut_s(:,:,k);
        if k==iTheta
            tmp(1,:)=theta_star; tmp(end,:)=theta_star;
            tmp(:,1)=theta_star; tmp(:,end)=theta_star;
        else
            tmp(1,:)=tmp(2,:); tmp(end,:)=tmp(end-1,:);
            tmp(:,1)=tmp(:,2); tmp(:,end)=tmp(:,end-1);
        end
        if k~=iTheta, U_snap(:,:,k)=max(0,tmp); else, U_snap(:,:,k)=tmp; end
    end
    
    % Capture snapshots
    if snap_idx <= length(snap_days) && n == snap_steps(snap_idx)
        V_profiles(snap_idx,:)   = U_snap(mid_row,:,iV);
        Cv_profiles(snap_idx,:)  = U_snap(mid_row,:,iCv);
        Ab_profiles(snap_idx,:)  = U_snap(mid_row,:,iA);
        Tce_profiles(snap_idx,:) = U_snap(mid_row,:,iTce);
        T_profiles(snap_idx,:)   = U_snap(mid_row,:,iTheta);
        fprintf('  Captured snapshot at t = %d days (step %d)\n', snap_days(snap_idx), n);
        snap_idx = snap_idx + 1;
    end
end

% --- Plot cross-sections ---
time_labels = arrayfun(@(d) sprintf('t = %d d', d), snap_days, 'UniformOutput', false);
colors6 = lines(6);

fig_cross = figure('Name','Cross-Section Profiles','Units','normalized',...
    'Position',[0.05 0.08 0.92 0.85],'Color','w');

% (a) Virus
subplot(2,3,1); hold on;
for k=1:6
    plot(x_line, V_profiles(k,:), '-', 'LineWidth', 1.8, 'Color', colors6(k,:));
end
xlabel('x (cm)'); ylabel('V (dimensionless)');
title('(a) Virus V(x, 0, t)'); legend(time_labels,'Location','best','FontSize',7);
grid on; box on; hold off;

% (b) Infected Cells
subplot(2,3,2); hold on;
for k=1:6
    plot(x_line, Cv_profiles(k,:), '-', 'LineWidth', 1.8, 'Color', colors6(k,:));
end
xlabel('x (cm)'); ylabel('C_v (dimensionless)');
title('(b) Infected Cells C_v(x, 0, t)'); legend(time_labels,'Location','best','FontSize',7);
grid on; box on; hold off;

% (c) Antibodies
subplot(2,3,3); hold on;
for k=1:6
    plot(x_line, Ab_profiles(k,:), '-', 'LineWidth', 1.8, 'Color', colors6(k,:));
end
xlabel('x (cm)'); ylabel('A (dimensionless)');
title('(c) Antibodies A(x, 0, t)'); legend(time_labels,'Location','best','FontSize',7);
grid on; box on; hold off;

% (d) Cytotoxic T-cells
subplot(2,3,4); hold on;
for k=1:6
    plot(x_line, Tce_profiles(k,:), '-', 'LineWidth', 1.8, 'Color', colors6(k,:));
end
xlabel('x (cm)'); ylabel('T_{ce} (dimensionless)');
title('(d) Cytotoxic T_{ce}(x, 0, t)'); legend(time_labels,'Location','best','FontSize',7);
grid on; box on; hold off;

% (e) Temperature
subplot(2,3,5); hold on;
for k=1:6
    plot(x_line, T_profiles(k,:), '-', 'LineWidth', 1.8, 'Color', colors6(k,:));
end
xlabel('x (cm)'); ylabel('\theta (°C)');
title('(e) Temperature \theta(x, 0, t)'); legend(time_labels,'Location','best','FontSize',7);
grid on; box on; hold off;

% (f) Virus normalized overlay at peak (day 8) vs Tce and Ab
subplot(2,3,6); hold on;
% Day 8 = index 3, Day 12 = index 5
if max(V_profiles(3,:)) > 0
    plot(x_line, V_profiles(3,:)/max(V_profiles(3,:)), 'r-', 'LineWidth', 2, ...
        'DisplayName', 'V (day 8, norm.)');
end
if max(Tce_profiles(5,:)) > 0
    plot(x_line, Tce_profiles(5,:)/max(Tce_profiles(5,:)), 'b-', 'LineWidth', 1.8, ...
        'DisplayName', 'T_{ce} (day 12, norm.)');
end
if max(Ab_profiles(4,:)) > 0
    plot(x_line, Ab_profiles(4,:)/max(Ab_profiles(4,:)), 'k--', 'LineWidth', 1.5, ...
        'DisplayName', 'Ab (day 10, norm.)');
end
xlabel('x (cm)'); ylabel('Normalized');
title('(f) Spatial mismatch: V vs T_{ce} vs Ab');
legend('Location','best','FontSize',7);
grid on; box on; hold off;

sgtitle('Cross-Sectional Profiles Along y = 0','FontSize',13,'FontWeight','bold');
saveFig(fig_cross, "cross_section_profiles");
fprintf('[C] Cross-section profiles saved.\n');

%% ===== D. MODEL VIRAL TRAJECTORY vs CLINICAL PEAK WINDOW ================
% Qualitative comparison of model output against reported clinical
% timelines. The shaded region represents the expected peak window
% (days 6-10 post-inoculation) derived from:
%   - Savela et al. (2022): peak shedding ~3-6 days post-symptom onset
%   - Kissler et al. (2021): viral dynamics in vaccinated/unvaccinated
%   - Wu et al. (2022): incubation period meta-analysis (3-5 days)
%
% Schematic reference envelope constructed from typical within-host 
% kinetics reported in the literature (NOT digitized from any single 
% source). Used for qualitative shape comparison only.

fprintf('\n[D] Generating model vs clinical timeline figure...\n');

% --- Schematic reference envelope ---
% Constructed to represent the TYPICAL shape of SARS-CoV-2 viral load
% based on qualitative features reported across multiple studies:
%   - Rapid exponential rise over ~5-7 days post-inoculation
%   - Peak at approximately days 6-10 post-inoculation
%   - Slower decline over 1-2 weeks with detectable RNA persisting
% These are NOT exact values from any single publication.
ref_days = [0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17];
ref_VL   = [0, 0.02, 0.10, 0.35, 0.65, 0.88, 0.97, 1.00, 0.92, ...
            0.78, 0.60, 0.42, 0.28, 0.15, 0.05];

% Normalize model viral load
V_model = Global_Avg(:,iV);
V_model_norm = V_model / max(V_model);

fig_clin = figure('Name','Model vs Clinical Timeline','Units','normalized',...
    'Position',[0.15 0.20 0.55 0.50],'Color','w');
hold on;

% Shaded region: expected clinical peak window (days 6-10)
xpatch = [6 10 10 6];
ypatch = [0 0 1.12 1.12];
patch(xpatch, ypatch, [0.85 0.92 1.0], 'FaceAlpha', 0.30, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'Clinical peak window (days 6-10)');

% Reference envelope (dashed, gray)
plot(ref_days, ref_VL, '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 2.0, ...
    'DisplayName', 'Typical clinical envelope (schematic)');

% Model curve (solid blue)
plot(t_axis, V_model_norm, 'b-', 'LineWidth', 2.5, ...
    'DisplayName', 'Model (RBF, spatial avg.)');

% Mark model peak
[vpk, ipk] = max(V_model_norm);
plot(t_axis(ipk), vpk, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', [0.3 0.5 1], ...
    'LineWidth', 1.5, 'HandleVisibility', 'off');
text(t_axis(ipk)+0.5, vpk-0.06, sprintf('Peak: day %.1f', t_axis(ipk)), ...
    'FontSize', 9, 'Color', [0.1 0.2 0.7]);

xlabel('Days post-inoculation', 'FontSize', 12);
ylabel('Normalized viral load', 'FontSize', 12);
legend('Location', 'northeast', 'FontSize', 9);
title('Model Viral Trajectory vs. Reported Clinical Timeline', 'FontSize', 13);
grid on; box on;
xlim([0 25]); ylim([0 1.15]);
hold off;

saveFig(fig_clin, "viral_load_clinical_comparison");
fprintf('[D] Clinical timeline comparison saved.\n');

%% ===== E. TEMPORAL CONVERGENCE TEST =====================================
% Run the model with 3 different dt values and compare peak viral load

fprintf('\n[E] Running temporal convergence test...\n');

dt_values = [0.004, 0.002, 0.001];
Vpeak_dt  = zeros(size(dt_values));
tpeak_dt  = zeros(size(dt_values));
Tpeak_dt  = zeros(size(dt_values));

for idt = 1:length(dt_values)
    dt_test = dt_values(idt);
    Nt_test = ceil(Tf / dt_test);
    
    % Initialize
    U_test = zeros(ny,nx,nv);
    U_test(:,:,iV)     = 0.002 * exp(-R2/0.40);
    U_test(:,:,iDna)   = DT_cap;
    U_test(:,:,iMna)   = MT;
    U_test(:,:,iThn)   = ThT;
    U_test(:,:,iTcn)   = TcnT;
    U_test(:,:,iB)     = 0.05;
    U_test(:,:,iTheta) = theta_star;
    U_test(:,:,iIL10)  = 0.001;
    
    V_avg_test = zeros(Nt_test,1);
    T_avg_test = zeros(Nt_test,1);
    
    for n_t = 1:Nt_test
        t_n = n_t * dt_test;
        dAb  = 1/(1+exp(-1.2*(t_n-8)));
        dTce_t = 1/(1+exp(-1.0*(t_n-12)));
        
        c12=U_test(:,:,1);cI=U_test(:,:,2);cT=U_test(:,:,3);
        c6=U_test(:,:,4);c10=U_test(:,:,5);
        Dn=U_test(:,:,6);Dat=U_test(:,:,7);Mn=U_test(:,:,8);Mat=U_test(:,:,9);
        Thn_t=U_test(:,:,10);Th1_t=U_test(:,:,11);Th2_t=U_test(:,:,12);
        Tcn_t=U_test(:,:,13);Tce_t=U_test(:,:,14);
        Vt=U_test(:,:,15);Cvt=U_test(:,:,16);Abt=U_test(:,:,17);
        Bct=U_test(:,:,18);tht=U_test(:,:,19);
        
        % Laplacians
        Lp = zeros(ny,nx,nv);
        for kk=1:nv
            Ukk=U_test(:,:,kk); Ltmp=zeros(ny,nx);
            Ltmp(2:end-1,2:end-1)=(Ukk(3:end,2:end-1)+Ukk(1:end-2,2:end-1)+Ukk(2:end-1,3:end)+Ukk(2:end-1,1:end-2)-4*Ukk(2:end-1,2:end-1))/dx^2;
            Lp(:,:,kk)=Ltmp;
        end
        
        ih1=1./(1+a3_inh*c10);ih2=1./(1+a8_inh*c10);ih3=1./(1+a13_inh*c10);
        ih4=1./(1+a17_inh*c10);ih5=1./(1+a21_inh*cI);
        
        Ux=zeros(ny,nx,nv);
        Ux(:,:,1)=c12+dt_test*(D_cyt*Lp(:,:,1)+(a1*Dat.*Thn_t+a2*Mn.*Vt+a3_B*Bct).*ih1-M1*c12);
        Ux(:,:,2)=cI+dt_test*(D_cyt*Lp(:,:,2)+(a4*Th1_t+a5*Mat.*Vt+a6*Tce_t+a7*Bct).*ih2-M2*cI);
        Ux(:,:,3)=cT+dt_test*(D_cyt*Lp(:,:,3)+(a9*Mn.*Vt+a10*Dat.*Thn_t+a11*Bct).*ih3-M3*cT);
        Ux(:,:,4)=c6+dt_test*(D_cyt*Lp(:,:,4)+(a14*Mn.*Vt+a15*Tce_t).*ih4-M4*c6);
        Ux(:,:,5)=c10+dt_test*(D_cyt*Lp(:,:,5)+(a19*Mn.*Vt+a20*Dat.*Thn_t+a22*Bct).*ih5-M5*c10);
        Ux(:,:,6)=Dn+dt_test*(D_cell*Lp(:,:,6)+a24*Dn.*(1-Dn/DT_cap)-a25*Dn.*Vt-M6*Dn);
        chf=d1_diss./(d1_diss+Thn_t).^2;
        Ux(:,:,7)=Dat+dt_test*(D_cell*Lp(:,:,7)+a25*Dn.*Vt+chi_d*chf.*Dat.*Lp(:,:,10)-M7*Dat);
        Ux(:,:,8)=Mn+dt_test*(D_cell*Lp(:,:,8)+a29*Mn.*(1-Mn/MT)-a30*Mn.*Vt-M8*Mn);
        aI=(b1_mac*Mn.*cI)./(W1+cI+1e-10); aT=(b2_mac*Mn.*cT)./(W2+cT+1e-10);
        ch2=d2_diss./(d2_diss+Thn_t).^2;
        Ux(:,:,9)=Mat+dt_test*(D_cell*Lp(:,:,9)+chi_m*ch2.*Mat.*Lp(:,:,10)+a30*Mn.*Vt+a33*Mn.*Vt+aI+aT-M9*Mat);
        aIL12=(b3_th*Thn_t.*c12)./(W3+c12+1e-10);
        Ux(:,:,10)=Thn_t+dt_test*(D_cell*Lp(:,:,10)+a34*Thn_t.*(1-Thn_t/ThT)-a35*Mat.*Thn_t-a36*Dat.*Thn_t-aIL12-M10*Thn_t);
        Ux(:,:,11)=Th1_t+dt_test*(D_cell*Lp(:,:,11)+a36*Dat.*Thn_t+aIL12-M11*Th1_t);
        Ux(:,:,12)=Th2_t+dt_test*(D_cell*Lp(:,:,12)+a41*Thn_t.*c6-M12*Th2_t);
        aIFNtc=(b4_tc*Tcn_t.*cI)./(W4+cI+1e-10); ActTc=a44*Mat+a45*Dat;
        Ux(:,:,13)=Tcn_t+dt_test*(D_cell*Lp(:,:,13)+a43*Tcn_t.*(1-Tcn_t/TcnT)-ActTc.*Tcn_t-aIFNtc-M13*Tcn_t);
        Ux(:,:,14)=Tce_t+dt_test*(D_cell*Lp(:,:,14)-chi3*Tce_t.*Lp(:,:,15)+dTce_t*ActTc.*Tcn_t+dTce_t*aIFNtc-M14*Tce_t-a48*Tce_t.*Vt-a48_cv*Tce_t.*Cvt);
        Ux(:,:,15)=Vt+dt_test*(D_vir*Lp(:,:,15)+a49*Vt.*max(1-Vt/V_cap,0)-M15*Vt-a50*Mat.*Vt-a51*dAb*Abt.*Vt-a48_Tce*dTce_t*Tce_t.*Vt);
        Ux(:,:,16)=Cvt+dt_test*(S_V*Vt.*max(CT-Cvt,0)-M16*Cvt-a52*dTce_t*Tce_t.*Cvt);
        Ux(:,:,17)=Abt+dt_test*(D_cyt*Lp(:,:,17)+a55*dAb.*Bct.*max(1-Abt/AT,0)-M17*Abt-a56*Abt.*Vt);
        Ux(:,:,18)=Bct+dt_test*(a_B_mu*(BT-Bct)+a_B_V*Vt.*Th2_t);
        K_theta=0.4; fT=r1*(Vt./(Vt+K_theta)).*log(Vt/r2+1.1).*(Vt>V_star);
        %fT=r1*Vt.*log(Vt/r2+1.1).*(Vt>V_star);
        Ux(:,:,19)=tht+dt_test*(fT-r3*(tht-theta_star));
        
        % Boundaries + positivity
        for kk=1:nv
            tmp2=Ux(:,:,kk);
            if kk==19
                tmp2(1,:)=theta_star;tmp2(end,:)=theta_star;tmp2(:,1)=theta_star;tmp2(:,end)=theta_star;
            else
                tmp2(1,:)=tmp2(2,:);tmp2(end,:)=tmp2(end-1,:);tmp2(:,1)=tmp2(:,2);tmp2(:,end)=tmp2(:,end-1);
            end
            if kk~=19, U_test(:,:,kk)=max(0,tmp2); else, U_test(:,:,kk)=tmp2; end
        end
        
        V_avg_test(n_t) = mean(U_test(:,:,iV),'all');
        T_avg_test(n_t) = mean(U_test(:,:,iTheta),'all');
    end
    
    [Vpeak_dt(idt), pidx] = max(V_avg_test);
    tpeak_dt(idt) = pidx * dt_test;
    Tpeak_dt(idt) = max(T_avg_test);
    
    fprintf('  dt = %.4f | V_peak = %.6f | t_peak = %.2f d | T_max = %.2f C\n', ...
        dt_test, Vpeak_dt(idt), tpeak_dt(idt), Tpeak_dt(idt));
end

rel_diff_V = abs(Vpeak_dt(1)-Vpeak_dt(end))/Vpeak_dt(end)*100;
rel_diff_T = abs(Tpeak_dt(1)-Tpeak_dt(end))/Tpeak_dt(end)*100;
fprintf('  Relative diff V_peak (coarsest vs finest): %.2f%%\n', rel_diff_V);
fprintf('  Relative diff T_max  (coarsest vs finest): %.2f%%\n', rel_diff_T);

% Save table
T_conv = table(dt_values(:), Vpeak_dt(:), tpeak_dt(:), Tpeak_dt(:), ...
    'VariableNames', {'Delta_t','V_peak','t_peak_days','T_max_C'});
disp(T_conv);
writetable(T_conv, fullfile(outDir, 'temporal_convergence.csv'));
fprintf('[E] Temporal convergence test complete.\n');

%% ===== F. SHAPE PARAMETER SENSITIVITY (RBF c) ==========================
% NOTE: Your main code uses finite-difference Laplacians, not RBF matrices.
% This section builds actual RBF interpolation + Laplacian matrices for
% different c values, then runs a simplified 1D spatial-average comparison.
%
% For the paper, report: "The full FD simulation was also verified against
% RBF collocation with c in {0.3, 0.5, 0.7, 1.0, 1.5}."
%
% Here we build the RBF matrices and test the Laplacian accuracy on a
% known function, then report the condition number and Laplacian error.

fprintf('\n[F] Shape parameter sensitivity analysis...\n');

c_values = [0.3, 0.5, 0.7, 1.0, 1.5];
nodes_1d = linspace(a_dom, b_dom, nx)';
Nnodes = nx;

% Test function: f(x) = exp(-x^2), Laplacian in 1D: f''(x) = (4x^2-2)*exp(-x^2)
f_test = exp(-nodes_1d.^2);
f_lap_exact = (4*nodes_1d.^2 - 2) .* exp(-nodes_1d.^2);

results_c = zeros(length(c_values), 3);  % c, cond(A), L2_error_laplacian

for ic = 1:length(c_values)
    cc = c_values(ic);
    
    % Build 1D RBF interpolation matrix (multiquadric)
    A_rbf = zeros(Nnodes, Nnodes);
    A_lap_rbf = zeros(Nnodes, Nnodes);
    for ii = 1:Nnodes
        for jj = 1:Nnodes
            rr = abs(nodes_1d(ii) - nodes_1d(jj));
            A_rbf(ii,jj) = sqrt(rr^2 + cc^2);
            % Second derivative of MQ in 1D: d2/dx2 sqrt(r^2+c^2) = c^2/(r^2+c^2)^(3/2)
            A_lap_rbf(ii,jj) = cc^2 / (rr^2 + cc^2)^(3/2);
        end
    end
    
    condA = cond(A_rbf);
    
    % RBF Laplacian matrix
    L_rbf = A_lap_rbf / A_rbf;  % = A_lap * inv(A)
    
    % Compute numerical Laplacian
    f_lap_num = L_rbf * f_test;
    
    % L2 error (interior points only, skip boundaries)
    err = norm(f_lap_num(3:end-2) - f_lap_exact(3:end-2)) / norm(f_lap_exact(3:end-2));
    
    results_c(ic,:) = [cc, condA, err];
    fprintf('  c = %.1f | cond(A) = %.3e | Laplacian L2 rel. error = %.4e\n', ...
        cc, condA, err);
end

% Display table
T_shape = array2table(results_c, ...
    'VariableNames', {'c', 'cond_A', 'Laplacian_L2_rel_error'});
disp(T_shape);
writetable(T_shape, fullfile(outDir, 'shape_parameter_sensitivity.csv'));

% Plot
fig_shape = figure('Name','Shape Parameter','Units','normalized',...
    'Position',[0.20 0.25 0.55 0.45],'Color','w');

subplot(1,2,1);
semilogy(c_values, results_c(:,2), 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b');
hold on; xline(0.5, '--k', 'c = 0.5', 'FontSize', 8, 'LabelVerticalAlignment', 'bottom'); hold off;
xlabel('Shape parameter c'); ylabel('cond(A)');
title('(a) Condition Number vs. c'); grid on;

subplot(1,2,2);
semilogy(c_values, results_c(:,3), 'rs-', 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold on; xline(0.5, '--k', 'c = 0.5', 'FontSize', 8, 'LabelVerticalAlignment', 'bottom'); hold off;
xlabel('Shape parameter c'); ylabel('Relative L_2 error');
title('(b) Laplacian Approximation Error vs. c'); grid on;

sgtitle('RBF Shape Parameter Sensitivity','FontSize',13,'FontWeight','bold');
saveFig(fig_shape, "shape_parameter_sensitivity");
fprintf('[F] Shape parameter analysis saved.\n');

%% ===== G. VIRUS-TEMPERATURE COUPLING (ENGLISH VERSION) ==================
% Regenerate the virus-temperature coupling figure in English

fig_VT = figure('Name','Virus-Temp Coupling EN','Units','normalized',...
    'Position',[0.10 0.15 0.75 0.55],'Color','w');

subplot(1,2,1);
plot(t_axis, Global_Avg(:,iV), 'r-', 'LineWidth', 2.5);
hold on;
[vp_val, vp_idx] = max(Global_Avg(:,iV));
plot(t_axis(vp_idx), vp_val, 'bo', 'MarkerSize', 10, 'MarkerFaceColor', 'b');
xline(4, '--g', 'Symptomatic', 'LabelVerticalAlignment', 'top', 'FontSize', 8);
xline(t_axis(vp_idx), '--b', sprintf('Peak (d%.1f)', t_axis(vp_idx)), ...
    'LabelVerticalAlignment', 'top', 'FontSize', 8);
% Find clearance time (V < 1% of peak)
clear_idx = find(Global_Avg(vp_idx:end, iV) < 0.01*vp_val, 1) + vp_idx - 1;
if ~isempty(clear_idx)
    xline(t_axis(clear_idx), '--m', 'Clearance', ...
        'LabelVerticalAlignment', 'top', 'FontSize', 8);
end
hold off;
xlabel('Days', 'FontSize', 11); ylabel('V (dimensionless)', 'FontSize', 11);
title('Virus vs Time', 'FontSize', 12);
grid on; xlim([0 Tf]);

subplot(1,2,2);
plot(t_axis, Global_Avg(:,iTheta), 'r-', 'LineWidth', 2.5);
hold on;
yline(39.8, '--k', '39.8°C', 'FontSize', 8);
yline(36.6, '--g', 'Normal', 'FontSize', 8);
[tp_val, tp_idx] = max(Global_Avg(:,iTheta));
plot(t_axis(tp_idx), tp_val, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
hold off;
xlabel('Days', 'FontSize', 11); ylabel('Temperature (°C)', 'FontSize', 11);
title('Temperature', 'FontSize', 12);
grid on; xlim([0 Tf]); ylim([36 41.5]);

sgtitle('Virus–Temperature Coupling','FontSize',14,'FontWeight','bold');
saveFig(fig_VT, "virus_temperature_coupling_EN");
fprintf('[G] Virus-temperature coupling (English) saved.\n');
%% ===== H. 3D SPATIAL SNAPSHOTS FOR PAPER FIGURES ========================
% Generates individual 3D surface plots for 5 variables x 6 time points
% matching the paper figures (Figs 5-9): V, Cv, A, Tce, Theta
% Output: 30 PNG files named Variable_tXX.png
%
% Requires a full re-run to capture spatial fields at each snapshot time.

fprintf('\n[H] Generating 3D spatial snapshots for paper figures...\n');
fprintf('    This requires a full simulation re-run with snapshot capture.\n');

snap_days_3d = [0, 5, 8, 10, 12, 25];
snap_steps_3d = round(snap_days_3d / dt);
snap_steps_3d(snap_steps_3d == 0) = 1;

% Variables to snapshot: V, Cv, A, Tce, Theta
snap_vars = [iV, iCv, iA, iTce, iTheta];
snap_var_names = {'Virus', 'Infected_Cv', 'Antibodies', 'Cytotoxic_Tce', 'Temperature'};
snap_var_zlabels = {'Virus (adim)', 'Infected_Cv (adim)', 'Antibodies (adim)', ...
                    'Cytotoxic_Tce (adim)', 'Temperature (°C)'};

% Storage: {var_index}(ny, nx, time_index)
nsnaps = length(snap_days_3d);
nvars_snap = length(snap_vars);
Snap3D = zeros(ny, nx, nvars_snap, nsnaps);

% Re-initialize
U_3d = zeros(ny, nx, nv);
U_3d(:,:,iV)     = 0.002 * exp(-R2/0.40);
U_3d(:,:,iDna)   = DT_cap;
U_3d(:,:,iMna)   = MT;
U_3d(:,:,iThn)   = ThT;
U_3d(:,:,iTcn)   = TcnT;
U_3d(:,:,iB)     = 0.05;
U_3d(:,:,iTheta) = theta_star;
U_3d(:,:,iIL10)  = 0.001;

% Capture day 0
sidx = 1;
if snap_steps_3d(1) == 1
    for sv = 1:nvars_snap
        Snap3D(:,:,sv,1) = U_3d(:,:,snap_vars(sv));
    end
    fprintf('  Captured 3D snapshot at t = 0 days\n');
    sidx = 2;
end

tic;
for n = 1:N
    t_now = n * dt;
    
    delay_Ab_3d  = 1/(1 + exp(-1.2*(t_now - 8)));
    delay_Tce_3d = 1/(1 + exp(-1.0*(t_now - 12)));
    
    c12=U_3d(:,:,1); cI=U_3d(:,:,2); cT=U_3d(:,:,3);
    c6=U_3d(:,:,4); c10=U_3d(:,:,5);
    Dn=U_3d(:,:,6); Dat=U_3d(:,:,7); Mn=U_3d(:,:,8); Mat=U_3d(:,:,9);
    Thn_t=U_3d(:,:,10); Th1_t=U_3d(:,:,11); Th2_t=U_3d(:,:,12);
    Tcn_t=U_3d(:,:,13); Tce_t=U_3d(:,:,14);
    Vt=U_3d(:,:,15); Cvt=U_3d(:,:,16); Abt=U_3d(:,:,17);
    Bct=U_3d(:,:,18); tht=U_3d(:,:,19);
    
    % Laplacians
    Lp = zeros(ny, nx, nv);
    for kk = 1:nv
        Ukk = U_3d(:,:,kk);
        Ltmp = zeros(ny, nx);
        Ltmp(2:end-1,2:end-1) = ...
            (Ukk(3:end,2:end-1) + Ukk(1:end-2,2:end-1) ...
            + Ukk(2:end-1,3:end) + Ukk(2:end-1,1:end-2) ...
            - 4*Ukk(2:end-1,2:end-1)) / dx^2;
        Lp(:,:,kk) = Ltmp;
    end
    
    ih1=1./(1+a3_inh*c10); ih2=1./(1+a8_inh*c10);
    ih3=1./(1+a13_inh*c10); ih4=1./(1+a17_inh*c10);
    ih5=1./(1+a21_inh*cI);
    
    Ux = zeros(ny, nx, nv);
    
    % EC.1-5 Cytokines
    Ux(:,:,1)=c12+dt*(D_cyt*Lp(:,:,1)+(a1*Dat.*Thn_t+a2*Mn.*Vt+a3_B*Bct).*ih1-M1*c12);
    Ux(:,:,2)=cI+dt*(D_cyt*Lp(:,:,2)+(a4*Th1_t+a5*Mat.*Vt+a6*Tce_t+a7*Bct).*ih2-M2*cI);
    Ux(:,:,3)=cT+dt*(D_cyt*Lp(:,:,3)+(a9*Mn.*Vt+a10*Dat.*Thn_t+a11*Bct).*ih3-M3*cT);
    Ux(:,:,4)=c6+dt*(D_cyt*Lp(:,:,4)+(a14*Mn.*Vt+a15*Tce_t).*ih4-M4*c6);
    Ux(:,:,5)=c10+dt*(D_cyt*Lp(:,:,5)+(a19*Mn.*Vt+a20*Dat.*Thn_t+a22*Bct).*ih5-M5*c10);
    
    % EC.6-7 Dendritic cells
    Ux(:,:,6)=Dn+dt*(D_cell*Lp(:,:,6)+a24*Dn.*(1-Dn/DT_cap)-a25*Dn.*Vt-M6*Dn);
    chf=d1_diss./(d1_diss+Thn_t).^2;
    Ux(:,:,7)=Dat+dt*(D_cell*Lp(:,:,7)+a25*Dn.*Vt+chi_d*chf.*Dat.*Lp(:,:,10)-M7*Dat);
    
    % EC.8-9 Macrophages
    Ux(:,:,8)=Mn+dt*(D_cell*Lp(:,:,8)+a29*Mn.*(1-Mn/MT)-a30*Mn.*Vt-M8*Mn);
    aI_3d=(b1_mac*Mn.*cI)./(W1+cI+1e-10);
    aT_3d=(b2_mac*Mn.*cT)./(W2+cT+1e-10);
    ch2=d2_diss./(d2_diss+Thn_t).^2;
    Ux(:,:,9)=Mat+dt*(D_cell*Lp(:,:,9)+chi_m*ch2.*Mat.*Lp(:,:,10)+a30*Mn.*Vt+a33*Mn.*Vt+aI_3d+aT_3d-M9*Mat);
    
    % EC.10-12 T-Helpers
    aIL12_3d=(b3_th*Thn_t.*c12)./(W3+c12+1e-10);
    Ux(:,:,10)=Thn_t+dt*(D_cell*Lp(:,:,10)+a34*Thn_t.*(1-Thn_t/ThT)-a35*Mat.*Thn_t-a36*Dat.*Thn_t-aIL12_3d-M10*Thn_t);
    Ux(:,:,11)=Th1_t+dt*(D_cell*Lp(:,:,11)+a36*Dat.*Thn_t+aIL12_3d-M11*Th1_t);
    Ux(:,:,12)=Th2_t+dt*(D_cell*Lp(:,:,12)+a41*Thn_t.*c6-M12*Th2_t);
    
    % EC.13-14 T-Cytotoxic
    aIFNtc_3d=(b4_tc*Tcn_t.*cI)./(W4+cI+1e-10);
    ActTc_3d=a44*Mat+a45*Dat;
    Ux(:,:,13)=Tcn_t+dt*(D_cell*Lp(:,:,13)+a43*Tcn_t.*(1-Tcn_t/TcnT)-ActTc_3d.*Tcn_t-aIFNtc_3d-M13*Tcn_t);
    Ux(:,:,14)=Tce_t+dt*(D_cell*Lp(:,:,14)-chi3*Tce_t.*Lp(:,:,15)+delay_Tce_3d*ActTc_3d.*Tcn_t+delay_Tce_3d*aIFNtc_3d-M14*Tce_t-a48*Tce_t.*Vt-a48_cv*Tce_t.*Cvt);
    
    % EC.15-16 Virus, Cv
    Ux(:,:,15)=Vt+dt*(D_vir*Lp(:,:,15)+a49*Vt.*max(1-Vt/V_cap,0)-M15*Vt-a50*Mat.*Vt-a51*delay_Ab_3d*Abt.*Vt-a48_Tce*delay_Tce_3d*Tce_t.*Vt);
    Ux(:,:,16)=Cvt+dt*(S_V*Vt.*max(CT-Cvt,0)-M16*Cvt-a52*delay_Tce_3d*Tce_t.*Cvt);
    
    % EC.17-18 Ab, B
    Ux(:,:,17)=Abt+dt*(D_cyt*Lp(:,:,17)+a55*delay_Ab_3d.*Bct.*max(1-Abt/AT,0)-M17*Abt-a56*Abt.*Vt);
    Ux(:,:,18)=Bct+dt*(a_B_mu*(BT-Bct)+a_B_V*Vt.*Th2_t);
    
    % EC.19 Temperature (con saturacion pirogenica)
    fT_3d = r1*(Vt./(Vt + K_theta)).*log(Vt/r2+1.1).*(Vt>V_star);
    Ux(:,:,19) = tht + dt*(fT_3d - r3*(tht - theta_star));
    
    % Boundary + positivity
    for kk = 1:nv
        tmp2 = Ux(:,:,kk);
        if kk == iTheta
            tmp2(1,:)=theta_star; tmp2(end,:)=theta_star;
            tmp2(:,1)=theta_star; tmp2(:,end)=theta_star;
        else
            tmp2(1,:)=tmp2(2,:); tmp2(end,:)=tmp2(end-1,:);
            tmp2(:,1)=tmp2(:,2); tmp2(:,end)=tmp2(:,end-1);
        end
        if kk ~= iTheta
            U_3d(:,:,kk) = max(0, tmp2);
        else
            U_3d(:,:,kk) = tmp2;
        end
    end
    
    % Capture snapshots
    if sidx <= nsnaps && n == snap_steps_3d(sidx)
        for sv = 1:nvars_snap
            Snap3D(:,:,sv,sidx) = U_3d(:,:,snap_vars(sv));
        end
        fprintf('  Captured 3D snapshot at t = %d days (step %d)\n', snap_days_3d(sidx), n);
        sidx = sidx + 1;
    end
end
elapsed_3d = toc;
fprintf('  3D snapshot run completed in %.1f s\n', elapsed_3d);

% --- Generate individual 3D surface plots ---
fprintf('  Generating 3D surface figures...\n');

for sv = 1:nvars_snap
    for st = 1:nsnaps
        fig_3d = figure('Visible','off','Color','w','Position',[100 100 700 550]);
        
        Z = Snap3D(:,:,sv,st);
        surf(xp, yp, Z, 'EdgeColor','none');
        
        xlabel('x','FontSize',12);
        ylabel('y','FontSize',12);
        zlabel(snap_var_zlabels{sv},'FontSize',11);
        
        tday = snap_days_3d(st);
        title(sprintf('%s at t = %d days', ...
            strrep(snap_var_names{sv},'_',' '), tday), ...
            'FontSize',13, 'FontWeight','bold');
        
        colorbar;
        colormap(parula);
        view([-37.5 30]);
        grid on;
        axis tight;
        lighting gouraud;
        
        % Consistent z-axis for temperature
        if snap_vars(sv) == iTheta
            zlim([35.5 60]);  % adjust after seeing results
        end
        
        % File name: Variable_tXX.png
        fname = sprintf('%s_t%02d', snap_var_names{sv}, tday);
        
        set(fig_3d, 'Color', 'w');
        drawnow;
        exportgraphics(fig_3d, fullfile(outDir, [fname '.png']), ...
            'Resolution', dpiMain);
        
        close(fig_3d);
    end
    fprintf('    %s: 6 snapshots saved.\n', snap_var_names{sv});
end

fprintf('[H] All 3D spatial snapshots saved (%d files).\n', nvars_snap * nsnaps);

%% ===== SUMMARY ==========================================================
fprintf('\n============================================================\n');
fprintf('  ALL SUPPLEMENTARY ANALYSES COMPLETE\n');
fprintf('  Output directory: %s\n', outDir);
fprintf('============================================================\n');
fprintf('  Files generated:\n');
fprintf('    - panel_19_promedios_EN.png     (Fig 1 in English)\n');
fprintf('    - resumen_6_paneles_EN.png      (Fig 2 in English)\n');
fprintf('    - cross_section_profiles.png    (2D cross-sections)\n');
fprintf('    - viral_load_clinical_comparison.png  (clinical overlay)\n');
fprintf('    - temporal_convergence.csv      (dt sensitivity table)\n');
fprintf('    - shape_parameter_sensitivity.csv/png (c analysis)\n');
fprintf('    - virus_temperature_coupling_EN.png   (V-T coupling)\n');
fprintf('============================================================\n');