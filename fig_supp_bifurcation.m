
clear; 
clc;
close all



addpath( ...
    'dde_biftool_v3.1.1/ddebiftool', ...
    'dde_biftool_v3.1.1/ddebiftool_extra_psol', ...
    'dde_biftool_v3.1.1/ddebiftool_utilities');





cv_values = [1.0 1.5];
wcc_min = 1.5;
wcc_max = 2.5;
M = 20;

ind_wcc = 9;
ind_tau = 21;


a_e = 0.7;      
a_i = 0.7;
a_T = 0.07;     
a_r = 0.035;
I_e = -0.5;     
I_i = -0.61;
I_T = 0.16;     
I_r = -0.514;
wee = 0.9;      
wei = 1.5;
wie = -2.2;     
wii = -0.5;
wer = 2.0;      
weT = 1.9;
wTi = 1.8;      
wTe = 1.5;
wrT = -1.6;     
wTr = 1.44;

figure();
tl = tiledlayout(1,numel(cv_values),'TileSpacing','compact','Padding','compact');

for cv_idx = 1:numel(cv_values)

    cv = cv_values(cv_idx);

    tau_callosal = zeros(1,M);
    for j = 1:M
        tau_value = 25;
        while tau_value >= 25
            tau_value = 4.5/abs(cv + 0.2*randn);
        end
        tau_callosal(j) = 0.1*tau_value;
    end


    par0 = [a_e, a_i, a_T, a_r, ...
            I_e, I_i, I_T, I_r, ...
            wcc_min, wee, wei, wie, wii, ...
            wer, weT, wTi, wTe, wrT, wTr, ...
            M, 0.5, tau_callosal];

    ind_delays = ind_tau:(ind_tau+M);


    funcs = set_funcs( ...
        'sys_rhs', @sys_rhs, ...
        'sys_tau', @() ind_delays, ...
        'sys_deri', @sys_deri, ...
        'x_vectorized', true);

    %% Steady-state branch %%
    guess = [-.26 -.26 -.33 -.33 -.068 -.068 -.073 -.073]';

    stst_method = df_mthod(funcs,'stst',1).point;
    stst_method.newton_max_iterations = 20;

    stst1 = struct('kind','stst','parameter',par0,'x',guess);
    [stst1,~] = p_correc(funcs,stst1,[],[],stst_method);

    stst2 = stst1;
    stst2.parameter(ind_wcc) = stst2.parameter(ind_wcc) + 0.01;
    [stst2,~] = p_correc(funcs,stst2,[],[],stst_method);

    br = df_brnch(funcs,ind_wcc,'stst');
    br.point = [stst1 stst2];
    br.parameter.free = ind_wcc;
    br.parameter.min_bound(1,:) = [ind_wcc wcc_min];
    br.parameter.max_bound(1,:) = [ind_wcc wcc_max];
    br.parameter.max_step(1,:) = [ind_wcc 0.005];
    br.method.continuation.detect_bifurcations = 1;
    br.method.continuation.plot = 0;
    br.method.stability.minimal_real_part = -5;

    br = br_contn(funcs,br,600);
    [nunst_eq,~,~,br.point] = GetStability(br,'funcs',funcs);
    eq_un = nunst_eq(:).' > 0;


    %% Periodic orbit%%
    
    indx_hopf = find(strcmp({br.point.kind},'hopf'),1);

    if isempty(indx_hopf)
        indx_hopf = find(eq_un,1);
    end
    
    
    hopf_pt = p_tohopf(funcs,br.point(indx_hopf));
    unstable = find(eq_un);
    indx_sim = unstable(round(numel(unstable)/2));
    stst_sim = br.point(indx_sim);
    
    par = stst_sim.parameter;
    x0 = stst_sim.x;
    
    dde_rhs = @(~,y,Z) sys_rhs([y,Z],par);
    
    t_end = 3000;
    sol = dde23(dde_rhs,par(21:end),@(t)x0-0.001,[0 t_end], ddeset('RelTol',1e-5,'AbsTol',1e-7));
    
    late = sol.x > 0.75*t_end;
    t = sol.x(late);
    E = sol.y(1,late);
    
    peaks = find(E(2:end-1)>E(1:end-2) & ...
                 E(2:end-1)>=E(3:end)) + 1;
        
    assert(numel(peaks)>=2)
    
    t0 = t(peaks(end-1));
    T_guess = t(peaks(end))-t0;
    
    t_sample = linspace(0,T_guess,1000);
    prof_sim = deval(sol,t0+t_sample);
    
    psol_method = df_mthod(funcs,'psol').point;
    psol_method.newton_max_iterations = 200;
    

    degree = 3;
    intervals = 60;
    psol0 = p_topsol(funcs,hopf_pt,1e-3,degree,intervals);
    psol0.parameter = par;
    psol0.period = T_guess;
    psol0.profile = interp1(t_sample,prof_sim.', ...
                            psol0.mesh*T_guess,'spline').';
    
    [psol1,~] = p_correc(funcs,psol0,[],[],psol_method);
    
    psol2 = psol1;
    psol2.parameter(ind_wcc) = psol2.parameter(ind_wcc) + 1e-3;
    
    [psol2,~] = p_correc(funcs,psol2,[],[],psol_method);

    %% Periodic-orbit branch %%
    pbr = df_brnch(funcs,ind_wcc,'psol');
    pbr.point = [psol1 psol2];
    pbr.parameter.free = ind_wcc;
    pbr.parameter.min_bound(1,:) = [ind_wcc wcc_min];
    pbr.parameter.max_bound(1,:) = [ind_wcc wcc_max];
    pbr.parameter.max_step(1,:) = [ind_wcc 5e-3];
    pbr.method.point = psol_method;
    pbr.method.continuation.plot = 0;

    pbr = br_contn(funcs,pbr,250);
    pbr = br_rvers(pbr);
    pbr = br_contn(funcs,pbr,250);

    amp_all = arrayfun(@(pt) ...
    max(pt.profile(1,:))-min(pt.profile(1,:)),pbr.point);

    keep = amp_all > 1e-5;

    
    pbr.point = pbr.point(keep);


    %% Floquet multipliers %%
    [nunst_po,~,~,pbr.point] = GetStability( ...
        pbr,'funcs',funcs,'exclude_trivial',true);

    stable_po = nunst_po(:).' == 0;
    unstable_po = nunst_po(:).' > 0;
    unresolved_po = isnan(nunst_po(:).');


    %% Plot %%
    
    ax = nexttile(tl);
    hold(ax,'on')
    box(ax,'on')
    
    set(ax,'FontSize',11, ...
           'LineWidth',0.8, ...
           'TickLabelInterpreter','latex')
    
    
    wcc_eq = arrayfun(@(pt) pt.parameter(ind_wcc),br.point);
    E_eq = arrayfun(@(pt) pt.x(1),br.point);
    
    E_stable = E_eq;
    E_stable(eq_un) = NaN;
    
    E_unstable = E_eq;
    E_unstable(~eq_un) = NaN;
    
    plot(ax,wcc_eq,E_stable,'-', ...
        'Color','g', ...
        'LineWidth',1.6, ...
        'DisplayName','stable fixed point');
    
    plot(ax,wcc_eq,E_unstable,'-', ...
        'Color','r', ...
        'LineWidth',1.6, ...
        'DisplayName','unstable fixed point');
    
    
    wcc_po = arrayfun(@(pt) pt.parameter(ind_wcc),pbr.point);
    maxE = arrayfun(@(pt) max(pt.profile(1,:)),pbr.point);
    minE = arrayfun(@(pt) min(pt.profile(1,:)),pbr.point);
    

    
    maxE_stable = maxE;
    minE_stable = minE;
    maxE_stable(~stable_po) = NaN;
    minE_stable(~stable_po) = NaN;
    
    maxE_unstable = maxE;
    minE_unstable = minE;
    maxE_unstable(~unstable_po) = NaN;
    minE_unstable(~unstable_po) = NaN;
    
    plot(ax,wcc_po,maxE_stable,'-', ...
        'Color','k', ...
        'LineWidth',1.2, ...
        'DisplayName','stable limit cycle');
    
    plot(ax,wcc_po,minE_stable,'-', ...
        'Color','k', ...
        'LineWidth',1.2, ...
        'HandleVisibility','off');
    
    plot(ax,wcc_po,maxE_unstable,'-', ...
        'Color','b', ...
        'LineWidth',1.2, ...
        'DisplayName','unstable limit cycle');
    
    plot(ax,wcc_po,minE_unstable,'-', ...
        'Color','b', ...
        'LineWidth',1.2, ...
        'HandleVisibility','off');
    
    if any(unresolved_po)
    
        maxE_unresolved = maxE;
        minE_unresolved = minE;
    
        maxE_unresolved(~unresolved_po) = NaN;
        minE_unresolved(~unresolved_po) = NaN;
    
        plot(ax,wcc_po,maxE_unresolved,'-', ...
            'Color',[0.60 0.60 0.60], ...
            'LineWidth',1.2, ...
            'DisplayName','unresolved');
    
        plot(ax,wcc_po,minE_unresolved,'-', ...
            'Color',[0.60 0.60 0.60], ...
            'LineWidth',1.2, ...
            'HandleVisibility','off');
    end
    
    
    xlim(ax,[wcc_min wcc_max])
    ylim(ax,[-1 0])
    
    xlabel(ax,'$w_{\mathrm{callosal}}$','Interpreter','latex')
    
    if cv_idx == 1
        ylabel(ax,'$V^{*}_{\mathrm{PYR}}$','Interpreter','latex')
    end
    
    title(ax,['$\mathrm{CV} = ' num2str(cv,'%.1f') '\ \mathrm{m/s}$'], ...
          'Interpreter','latex')
    
    legend(ax,'show', ...
           'Location','best', ...
           'Box','on', ...
           'Interpreter','latex', ...
           'FontSize',8)
end




function f = sys_rhs(xx,p)

    S = @(u) 1./(1+exp(-25*u));
    cidx = 3:(p(20)+2);

    a_e = p(1);  a_i = p(2);  a_T = p(3);  a_r = p(4);
    I_e = p(5);  I_i = p(6);  I_T = p(7);  I_r = p(8);
    wcc = p(9);  wee = p(10); wei = p(11); wie = p(12); wii = p(13);
    wer = p(14); weT = p(15); wTi = p(16); wTe = p(17);
    wrT = p(18); wTr = p(19); M = p(20);

    E1 = xx(1,1,:); E2 = xx(2,1,:);
    I1 = xx(3,1,:); I2 = xx(4,1,:);
    T1 = xx(5,1,:); T2 = xx(6,1,:);
    R1 = xx(7,1,:); R2 = xx(8,1,:);

    E1d = xx(1,2,:); E2d = xx(2,2,:);
    T1d = xx(5,2,:); T2d = xx(6,2,:);

    cc1 = sum(S(xx(2,cidx,:)),2);
    cc2 = sum(S(xx(1,cidx,:)),2);

    f = [a_e.*(-E1 + wee.*S(E1) + wie.*S(I1) + wTe.*S(T1d) ...
                    + (wcc/M).*cc1 + I_e); ...
         a_e.*(-E2 + wee.*S(E2) + wie.*S(I2) + wTe.*S(T2d) ...
                    + (wcc/M).*cc2 + I_e); ...
         a_i.*(-I1 + wei.*S(E1) + wii.*S(I1) + wTi.*S(T1d) + I_i); ...
         a_i.*(-I2 + wei.*S(E2) + wii.*S(I2) + wTi.*S(T2d) + I_i); ...
         a_T.*(-T1 + weT.*S(E1d) + wrT.*S(R1) + I_T); ...
         a_T.*(-T2 + weT.*S(E2d) + wrT.*S(R2) + I_T); ...
         a_r.*(-R1 + wTr.*S(T1) + wer.*S(E1d) + I_r); ...
         a_r.*(-R2 + wTr.*S(T2) + wer.*S(E2d) + I_r)];
end

function J = sys_deri(xx,p,nx,np,v)

if numel(nx) ~= 1 || ~isempty(np) || ~isempty(v)
    J = df_deriv(struct('sys_rhs',@sys_rhs),xx,p,nx,np,v);
    return
end

J = zeros(8,8,size(xx,3));

a_e = p(1);
a_i = p(2);
a_T = p(3);
a_r = p(4);

wcc = p(9);
wee = p(10);
wei = p(11);
wie = p(12);
wii = p(13);
wer = p(14);
weT = p(15);
wTi = p(16);
wTe = p(17);
wrT = p(18);
wTr = p(19);
M = p(20);

S = @(x) 1./(1+exp(-25*x));
dS = @(x) 25.*S(x).*(1-S(x));

if nx == 0

    E1 = xx(1,1,:);
    E2 = xx(2,1,:);
    I1 = xx(3,1,:);
    I2 = xx(4,1,:);
    T1 = xx(5,1,:);
    T2 = xx(6,1,:);
    R1 = xx(7,1,:);
    R2 = xx(8,1,:);

    J(1,1,:) = a_e.*(-1 + wee.*dS(E1));
    J(1,3,:) = a_e.*wie.*dS(I1);

    J(2,2,:) = a_e.*(-1 + wee.*dS(E2));
    J(2,4,:) = a_e.*wie.*dS(I2);

    J(3,1,:) = a_i.*wei.*dS(E1);
    J(3,3,:) = a_i.*(-1 + wii.*dS(I1));

    J(4,2,:) = a_i.*wei.*dS(E2);
    J(4,4,:) = a_i.*(-1 + wii.*dS(I2));

    J(5,5,:) = -a_T;
    J(5,7,:) = a_T.*wrT.*dS(R1);

    J(6,6,:) = -a_T;
    J(6,8,:) = a_T.*wrT.*dS(R2);

    J(7,5,:) = a_r.*wTr.*dS(T1);
    J(7,7,:) = -a_r;

    J(8,6,:) = a_r.*wTr.*dS(T2);
    J(8,8,:) = -a_r;

elseif nx == 1

    E1_delay = xx(1,2,:);
    E2_delay = xx(2,2,:);
    T1_delay = xx(5,2,:);
    T2_delay = xx(6,2,:);

    J(1,5,:) = a_e.*wTe.*dS(T1_delay);
    J(2,6,:) = a_e.*wTe.*dS(T2_delay);

    J(3,5,:) = a_i.*wTi.*dS(T1_delay);
    J(4,6,:) = a_i.*wTi.*dS(T2_delay);

    J(5,1,:) = a_T.*weT.*dS(E1_delay);
    J(6,2,:) = a_T.*weT.*dS(E2_delay);

    J(7,1,:) = a_r.*wer.*dS(E1_delay);
    J(8,2,:) = a_r.*wer.*dS(E2_delay);

else

    J(1,2,:) = (a_e*wcc/M).*dS(xx(2,nx+1,:));
    J(2,1,:) = (a_e*wcc/M).*dS(xx(1,nx+1,:));

end
end