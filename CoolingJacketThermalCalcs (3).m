clear; clc; clf;
psi2Pa = 6894.76;
%% Colors
Fuschia     = [168   1  99]/255;    %Establishing Colors
Violet      = [ 63  19 108]/255;
Dark_Blue   = [  8  66 126]/255;
Light_Blue  = [ 14 120 197]/255;
Turqoise    = [ 19 153 160]/255;
Light_Green = [103 180  25]/255;
Dark_Green  = [ 10  81  57]/255;
Dark_Grey   = [ 66  76  88]/255;
Red         = [237  28  36]/255;
Black       = [  0   0   0]/255;
Neon_Green  = [  0 255   0]/255;
Orange      = [255 155   0]/255;
IDK         = [150 150 150]/255;
Colors      = [Fuschia; Violet; Dark_Blue; Light_Blue; Turqoise; Light_Green; Dark_Green; Dark_Grey; Red; Black; Neon_Green; Orange; IDK]; 
%% AUF WIEDERSEHEN OURO

%INCONEL 718
%GAS Cp, k, and mu are from air data
%disk assumption for cooling channels
tic
%% Engine
g = 9.81; %gravity
epsilon = 1E-9; %tolerance to stop iterating within an axial step
nu = .92; %adiabatic wall temp stagnation recovery factor, assumed
num_disks = 100; %how many dx
er = 9.5e-6; %m
Rcarbon = 1/(5.678263398*144*3600*1/800); %thermal resistance of carbon deposits, 800 is from huzel then theres a bunch of messy conversion

%gas properties from cea
%Chamber
gc = 1.2252;
Rc = 451.2102464;
Pc = 500*psi2Pa; %Pa
Tc = 2803.74; %K combustion temp
Lc = .15; %m

%Throat
gt = 1.2395;
Rt = 450.0135318;
Lt = .05; %m

%Exit
ge = 1.2666;
Rexit = 449.7217432;
Le = .062; %m

%Engine Material and Geometry
throat_r = 0.02866370302/2; %throat radius, m
r_c = 0.085/2; %radius of chamber in meters
throat_curve = 0.01348627227;

%memory allocation for recrystalization temp
temp_limit = zeros(1, num_disks);
temp_limit(:) = 1273.15; %K

%axial steps
x = 0:(Lt+Lc+Le)/(num_disks-1):(Lt+Lc+Le);
dx = (Lt+Lc+Le)/(num_disks-1);

%engine geometry (m)
chamber = zeros(1, round(num_disks*Lc/(Lt+Lc+Le)));
chamber(:) = (r_c/throat_r)^2;
%area ratio for isentropic relations
area_ratio = [6.89:-(6.89-1)/((num_disks-1)*Le/(Lt+Lc+Le)):1 1:((r_c/throat_r)^2-1)/((num_disks-1)*Lt/(Lt+Lc+Le)):(r_c/throat_r)^2 chamber];
%equivalent diameter of channels
channel_size = [linspace(.0033, .00165,(num_disks)*Le/(Lt+Lc+Le)) .00165 linspace(.00165, .0015, (num_disks)*Lt/(Lt+Lc+Le)) ones(1, round(num_disks*Lc/(Lt+Lc+Le)))*.0015];
%inside wall thickness
THK = [.00176:-(.00176-.001)/((num_disks-1)*Le/(Lt+Lc+Le)):.001 linspace(.001, .001, (num_disks*(Lt+Lc)/(Lt+Lc+Le)))]+channel_size./2;

dh = 2*channel_size; %hydraulic diameter of cooling channel (disk assumption)

%Kero
mdot_K = 0.5130783892; %mass flow rate, kg/s
T_in_K = 298.15; %K

%combustion
gas_density_o = Pc/(Rc*Tc);
R = 450000; %SPECIFIC gas constant, J/kgK
mdot_g = .821+.53; %mas flow rate, kg/s
cstar = 1717.5; %m/s

gavg = 1.25305;
gamma = [ge, gt, gc];

%% Calc

%memory allocation of properties for each axial step
[qh, qw, qc, R_hot, R_wall, R_cold, h_cold, h_hot, Qdot, sigma, Tns, Thot, Tgv, Twh, Twc, Tcv, Tco, rhog, vg, cpg, kg, mu_g, prandtl_g, mu_k, rhok, cpk, vk, k_K, cpcoolant, rhocoolant, mucoolant, prandtl_c, M, A, r, channel_area, k_engine] ...
    = deal(ones(1, num_disks));
count = 0;

%loop through all dx
for i = 1:num_disks
    %get temp from previous step
    if i==1
        Tco(i) = T_in_K;
    else
        Tco(i) = Tco(i-1);
    end
    
    %to find mach number, findMach function is employed. have to tell it
    %whether the flow is supersonic or not so it uses the correct solution
    %diverging
    if i<=ceil((num_disks)*Le/(Lt+Lc+Le))
        loc = 1;%do this since CEA gives gamma at 3 points
        M(i) = findMach(area_ratio(i), gamma(loc), 2);
        
    %converging section
    elseif i < ceil(num_disks*Le/(Lt+Lc+Le)+num_disks*Lt/(Lt+Lc+Le))
        loc = 2;
        M(i) = findMach(area_ratio(i), gamma(loc), 1);
        
    %chamber
    else 
        loc = 3;
        M(i) = findMach(area_ratio(i), gamma(loc), 1);
    end
    
    A(i) = (pi*throat_r^2)*area_ratio(i);%compressible area relation to throat
    r(i) = sqrt(A(i)/pi);%inner radius of chamber
    rhog(i) = gas_density_o*(1+(gavg-1)/2*M(i)^2)^(-1/(gavg-1)); %compressible density relation
    Tns(i) = Tc*(1+(gavg-1)/2*M(i)^2)^(-1); %Compressible gas temp without correctional factor
    Thot(i) = Tns(i)*nu;%factor of safety on combustion temp (account for inefficiencies)
    vg(i) = sqrt(gavg*R*Thot(i)).*M(i); %mach number definition
    
    channel_area(i) = pi*((channel_size(i)+THK(i)+r(i))^2-(THK(i)+r(i))^2);
    
    %initialize guesses. Take hot side wall temp to be halfway between cold
    %and combustion. Take cold side wall to be 400 less than that (bc why
    %not). finally, get coolant temp from previous step
    Twh_guess = (Thot(i)+Tco(i))/2;
    Twc_guess = Twh_guess-400;
    Tco_guess = Tco(i);
    Qprev = 0;
    Q = 1;
    
    %while not converged
    while abs(Q-Qprev) > epsilon
        %% fluid properties calc
        %mostly just a bunch of polyfits of experimental data. got from ANSYS
        count = count + 1;
        Qprev = Q;
        
        %log mean across vapor boundaries
        Tgv(i) = (Thot(i)-Twh_guess)/(log(Thot(i)/Twh_guess));
        Tcv(i) = (Twc_guess-Tco_guess)/(log(Twc_guess/Tco_guess));
        
        %engine - https://www.hightempmetals.com/techdata/hitempInconel625data.php#5 data
        k_engine(i) = 0.0202*((Twh_guess-Twc_guess)/(log(Twh_guess/Twc_guess))) + 3.8872;
        
        %gas side
        %AIR FOR THESE THREE https://www.engineersedge.com/physics/viscosity_of_air_dynamic_and_kinematic_14483.htm#:~:text=The%20viscosity%20of%20air%20depends,the%20kinematic%20viscosity%2015.7%20cSt.
        cpg(i) = 1.112335e-19*Tgv(i)^7-1.565553e-15*Tgv(i)^6+9.237533e-12*Tgv(i)^5-2.936679e-08*Tgv(i)^4+5.421615e-05*Tgv(i)^3-0.0581276*Tgv(i)^2+33.70605*Tgv(i)-7069.814; 
        kg(i) = 4.64083*10^(-12)*Tgv(i)^3-2.52198*10^(-8)*Tgv(i)^2+7.56238*10^(-5)*Tgv(i)+.0236353;
        mu_g(i) = 1.22508E-15*Tgv(i)^3 - 9.03933E-12*Tgv(i)^2 + 3.65724E-08*Tgv(i) + 1.95125E-05;
        
        gas_Tr = Twh_guess/Tc; %gas side temperature ratio across boundary layer, uses stagnation temperature?
        sigma(i) = 1/(((.5*gas_Tr*(1+(gamma(loc)-1)/2*M(i)^2)+.5)^.68)*(1+(gamma(loc)-1)/2*M(i)^2)^.12);
        
        %coolant (bulk and vapor boundary properties calc'd)
        mu_k(i) = -1.0237e-15*Tcv(i)^5+(2.774e-12)*Tcv(i)^4+(-2.9991e-09)*Tcv(i)^3+(1.6207e-06)*Tcv(i)^2-0.00043949*Tcv(i)+0.048224;
        mucoolant(i) = (-1.0237e-15)*Tco_guess^5+(2.774e-12)*Tco_guess^4+(-2.9991e-09)*Tco_guess^3+(1.6207e-06)*Tco_guess^2-0.00043949*Tco_guess+0.048224;
        rhok(i) = -0.0002*Tcv(i)^2 - 0.5782*Tcv(i) + 995.72;
        rhocoolant(i) = -0.0002*Tco_guess^2 - 0.5782*Tco_guess + 995.72;
        cpk(i) = (0.0072*Tcv(i)-0.3671)*1000;
        cpcoolant(i) = (0.0072*Tco_guess-0.3671)*1000;
        k_K(i) = (-2.8346e-17)*Tcv(i)^3+(4.3498e-14)*Tcv(i)^2-0.0001756*Tcv(i)+0.172104;

        
        %% heat transfer section
        %calcing h_hot
        prandtl_g(i) = cpg(i)*mu_g(i)/kg(i);
        h_hot(i) = 1/(1/((.026/((2*throat_r)^.2))*((mu_g(i)^.2)*cpg(i)/(prandtl_g(i)^.6))*((Pc/cstar)^.8)*((2*throat_r/throat_curve)^.1)*((1/area_ratio(i))^.9)*sigma(i))+Rcarbon); %huzel 4-13
        
        %calcing h_cold
        prandtl_c(i) = mu_k(i)*cpk(i)/k_K(i);
        h_cold(i) = 1000*(.029*cpk(i)*(mu_k(i)^(.2))/(prandtl_c(i)^(2/3))*((mdot_K)^(.8)/dh(i)^(0.2))*((Tco_guess/Twc_guess)^(.55))); %huzel 4-25
        
        %thermal resistances
        R_hot(i) = 1/(h_hot(i)*2*pi*dx*r(i));
        R_wall(i) = THK(i)/(k_engine(i)*(2*pi*dx*((r(i)+THK(i))-r(i)))/log((r(i)+THK(i))/r(i)));
        R_cold(i) = 1/(h_cold(i)*2*pi*dx*(r(i)+THK(i)));

        %1 over R
        U = 1/(R_hot(i)+R_wall(i)+R_cold(i));
        
        qh(i) = h_hot(i)*(Thot(i)-Twh_guess);
        qw(i) = k_engine(i)/THK(i)*(Twh_guess-Twc_guess);
        qc(i) = h_cold(i)*(Twc_guess-Tco(i));
        
        %update coolant guess temp using arithmetic average across the dx
        if i == 1
            Tco_guess = .5*(Q/(mdot_K*cpcoolant(i)))+T_in_K;
        else
            Tco_guess = .5*(Q/(mdot_K*cpcoolant(i)))+Tco(i);
        end
        
        %iterational correction
        Twc_guess = Tco_guess+(Thot(i)-Tco_guess)*R_cold(i)/(1/U);
        Twh_guess = Twc_guess+(Thot(i)-Tco_guess)*R_wall(i)/(1/U);
        
        %total heat transfer according to guesses
        Q = U*(Thot(i)-Tco_guess);
    end
    %save the solution
    Tco(i) = Tco_guess;
    Twh(i) = Twh_guess;
    Twc(i) = Twc_guess;
    Qdot(i) = Q;
end

%plot
hold on
line = plot(x, Thot, x, Twh, x, Twc, x, Tco, x, temp_limit);
for ii = 1:length(line)
    line(ii).LineWidth = 2;
    line(ii).LineStyle = '-';
    line(ii).Marker    = 'none';
    line(ii).Color     = Colors(ii,:);
end
line(end).Color = Red;
title("Cooling Jacket Temperatures");
L = legend("Combustion Temp", "Hot Side Wall Temp", "Cold Side Wall Temp", "Coolant Temp", "Inconel 718 Recrystallization Temp", 'Location', 'best', 'FontSize', 14);
set(gcf, 'Color', 'white');
grid on;
grid minor;
ylabel("Temperature (K)");
xlabel("Distance from Nozzle Exit (m)");

figure(2)
plot(x, qh, x, qw, x, qc);
title("Cooling Jacket Heat Fluxi");
z = legend("qh", "qw", "qc");
set(gcf, 'Color', 'white');
grid on;
grid minor;
ylabel("Heat Flux (W/m^2-K)");
xlabel("Distance from Nozzle Exit (m)");
figure(7)
plot(x, Tco);

%% Delta P

vk = mdot_K./(channel_size.^2*pi*10.*rhocoolant);
t = dx./vk;
Re = rhocoolant.*vk.*dx./mucoolant;

%Moody ff first guess (turbulent flow)
f = 1.14+2.*log10(channel_size./er).^(-2);
fnew = ones(1, num_disks);
ffcount = 0;

%iteratively solve for ff using Colebrook Eq
while sum(sum(fnew-f)) > epsilon
    ffcount = ffcount + 1;
    f = fnew;
    fnew = (-2.*log10(((er./channel_size)./3.7)+(2.51./(Re.*(f.^.5))))).^(-2);
end

%Darcy Weisbach Eq
dP = dx*vk.^2./(channel_size).*.5.*rhocoolant.*fnew;
dPtot = sum(dP);
dPtotPSI = dPtot/psi2Pa;



static_stress = Pc*((Thot./Tc).^(gavg/(gavg-1))).*(r*2.+THK)./(2*THK);
figure(8)
plot(x, static_stress);
p_ratio = .28;
E = 184*10^9;
Sut = 630*10^6;
Sy = 280*10^6;
c = r+THK;
b = r;
alpha = 16*10^-6;
sigma_outer = (Twh-Twc)*alpha*E./(2.*(1-p_ratio).*log(c./b)).*(1-2*c.^2./(c.^2-b.^2).*log(c./b));
sigma_t = abs(sigma_outer)+static_stress;

FoS_b = Sut/max(sigma_t);
FoS_y = Sy/max(sigma_t);

toc





