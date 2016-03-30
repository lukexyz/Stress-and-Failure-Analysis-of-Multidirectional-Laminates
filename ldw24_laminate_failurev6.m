% ENME467 Computational Assignment
% Predicition of Laminate Failure using Classical Laminate Theory (CLT)
% By Luke Woods
% v6 April 2012

clear; clc; format shortg; format compact;

disp('---------------------------------------------------------------------------------------------- ');
disp('                            ENME467 CLT Laminate Failure Prediction                            ');
disp('                                         by Luke Woods                                         ');
disp('---------------------------------------------------------------------------------------------- ');
disp('Laminate_Failure v6                                                                  (May 2012)');
disp('Matlab R2011b')
disp(' ');
disp('For best viewing, please maximise the command window.');
disp(' ');
disp(' ');

% 1.0 Constants

% Table 1: Basic Engineering Constants
%           Vf     E1      E2      G12    v12
prop = [1, 0.55, 132.7E9, 8.83E9, 4.76E9, 0.36;
        2, 0.50, 043.0E9, 8.90E9, 4.50E9, 0.27;
        3, 0.58, 131.0E9, 8.70E9, 5.00E9, 0.28;
        4, 0.60, 087.0E9, 5.50E9, 2.20E9, 0.34;
        5, 0.50, 201.0E9, 21.7E9, 5.40E9, 0.17];

v21 = (prop(:,4).*prop(:,6))./prop(:,3);
% Table 2: Ultimate strengths (MPa)
prop2 = [1, 2280E6, 57E6, 1440E6, 228E6, 71E6;
         2, 1280E6, 49E6,  690E6, 158E6, 69E6;
         3, 2060E6, 78E6, 1080E6, 196E6, 157E6;
         4, 1280E6, 30E6,  335E6, 158E6, 49E6;
         5, 1380E6, 56E6, 1600E6, 125E6, 62E6];

% 1.1 Q matrices for each material (base 0-deg orientation)
Q = zeros(5*3,3);     
for z = 1:1:5
Q((z*3-2):(z*3),1:3) = [prop(z,3)/(1-(prop(z,6)*v21(z))), v21(z)*prop(z,3)/(1-(prop(z,6)*v21(z))), 0;
                        prop(z,6)*prop(z,4)/(1-prop(z,6)*v21(z)), prop(z,4)/(1-(prop(z,6)*v21(z))), 0;
                        0, 0, prop(z,5)];
end
          
          
% 2.0 ========  Input Variables ==========

% 2.1 Layer and thickness inputs
str = sprintf('      ENME467 CLT Laminate Failure Prediction \n                           by Luke Woods             \n'); 
str2 = sprintf('\n                                                                                 1. How many layers does the laminate have? ');
str = strcat(str, str2);
prompt = {str,' 2. What is the laminate thickness? (single ply, mm)'};
defaultanswer = {'8','0.19'};
name = 'CLT - Input Variables';
numlines = 1;
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
answer = inputdlg(prompt,name,numlines,defaultanswer,options);

lays = str2double(getfield(answer,{1}));
t = str2double(getfield(answer,{2}));
t = t/1000; % mm to m
Tt = t*lays; 


% 2.3 Initializing matrices pre for-loop
c = 0; % Current ply
plymat = zeros(lays*3,1); % Entries correspond to ply material
plytheta = zeros(lays*3,1); % Entries correspond to ply orientation
Tmatrix = zeros(lays*3,3);
Us = zeros(lays*3,5);

for i = 1:3:(lays*3)
    % 2.4 Material and orientation input
    c = c + 1;
    
    % Creating Input Box
    str1 = sprintf('============= Material Selection =============\n      1 = Carbon-epoxy       4 = Kevlar149-epoxy\n      2 = S-Glass-epoxy      5 = Boron-epoxy\n      3 = Carbon-PEEK         6 = Straw and Mud\n=========================================');
    str2 = sprintf('\n                                                                                     For ply number %d, please enter:', c);
    Q1 = sprintf('\nThe material type (1-6) ');
    Q2 = sprintf('The orientation (-90 to 90) ');
    str = strcat(str1, str2, Q1);
    prompt = {str,Q2};
    switch c
        case 1
            defaultanswer = {'1','0'};
        case 2
            defaultanswer = {'1','45'};
        case 3
            defaultanswer = {'1','-45'};
        case 4
            defaultanswer = {'1','90'};
        case 5
            defaultanswer = {'1','90'};
        case 6
            defaultanswer = {'1','-45'};
        case 7
            defaultanswer = {'1','45'};
        case 8
            defaultanswer = {'1','0'};
    end
    name = 'CLT - Input Variables';
    numlines = 1;
    options.Resize = 'on';
    options.WindowStyle = 'normal';
    options.Interpreter = 'tex';
    answer = inputdlg(prompt,name,numlines,defaultanswer,options);

    % 2.5 Extracting Inputs
    plymat(i) = str2double(getfield(answer,{1}));
    plytheta(i) = str2double(getfield(answer,{2}));
    while (plymat(i) ~= 1 && plymat(i) ~= 2 && plymat(i) ~= 3 && plymat(i) ~= 4 && plymat(i) ~= 5)
        % plymat(i) = input('--> Invalid choice. Please enter material type (1-5): ');
        prompt={'Invalid choice. Please choose a suitable material (1-5):'};
        name='Don''t be silly.';
        answer=inputdlg(prompt,name);
        plymat(i) = str2double(getfield(answer,{1}));
    end
            
    
   
    % 2.6 Transformation of Q0 matrix to Qbar
    m = cos(plytheta(i)*-1*pi/180);  % Cos(theta)
    n = sin(plytheta(i)*-1*pi/180) ; % Sin(theta)
    
    T = [m^2 n^2 -2*m*n; n^2 m^2 2*m*n; m*n -m*n m^2-n^2]; 
    Tmatrix(i:i+2,1:3) = T;
    
    Qbar(i:i+2,1:3) = T\Q(plymat(i):plymat(i)+2,1:3)/(T.');
    
    Us(i,1:5) = prop2(plymat(i,1),1:5);
end

% Redundant Waitbar
h = waitbar(0,'COMPUTING ADVANCED HYBRID ALGORITHM...');
steps = 600;
for step = 1:steps
    waitbar(step / steps)
end
close(h) 

% Laminate input display table
input1 = zeros(lays,3);
for f = 1:1:lays
    input1(f,1) = f;
    input1(f,2) = plytheta(f*3-2);
    input1(f,3) = plymat(f*3-2);
end

disp('User input 1:');
disp(' ' );
disp('--------------------------');
disp('     Laminate Config      ');
disp('--------------------------');
disp('   Ply   Dir   Material   ');
disp(input1);
disp('--------------------------');

str = sprintf('--> Target: %g ply laminate', lays);
disp(str);
disp('--> CLT calculations')
disp(' ');

% 3.0 Calculation of Laminate Stiffness Matrix (ABD)
h = -(lays/2)*t:t:(lays/2)*t;


%  =====Previously======
% Tt = total laminate thickness
% t = ply thickness

% 3.1 --- A matrix ---
A = zeros(3,3);
for i = 1:3;
    for j = 1:3;
       for k = 0:1:lays-1
            A(i,j) = A(i,j) + (Qbar(i+3*k,j)*(h(k+1)-h(k+2)));
       end
    end
end
A = (-1.*A);

% 3.2 --- B matrix ---
B = zeros(3,3);
% % % Check if symmetrical, 0 = symmetrical
% % sym = 0; 
% % for x = 1:1:lays/2
% %     if (input1(x,2) ~= input1(lays-x+1,2)) 
% %         sym = 1;
% %     end
% %     if (input1(x,3) ~= input1(lays-x+1,3))
% %         sym = 1;
% %     end
% % end


    for i = 1:3;
        for j = 1:3;
           for k = 0:1:lays-1
                B(i,j) = B(i,j) + (Qbar(i+3*k,j)*((h(k+1)^2)-(h(k+2)^2)));
           end
        end
    end
    B = (-1.*B/2);


% 3.3 --- D matrix ---
D = zeros(3,3);
for i = 1:3;
    for j = 1:3;
       for k = 0:1:lays-1
            D(i,j) = D(i,j) + (Qbar(i+3*k,j)*((h(k+1)^3)-(h(k+2)^3)));
       end
    end
end
D = (-1.*D/3);

% 3.4 Assembling and displaying ABD matrix
ABD = [A,B;B,D];
disp('================================ Laminate Stiffness Matrix (ABD) =============================');
disp(ABD);
disp('==============================================================================================');
disp(' ');

% 4.0 Forming Inverted Laminate Stiffness Matrix (abcd)

Bstar = -inv(A)*B;
Cstar = B/A;
Dstar = D - (B*inv(A))*B;
a = inv(A) - (Bstar*inv(Dstar))*Cstar;
b = Bstar*inv(Dstar);
c = -(inv(Dstar)*Cstar);
d = inv(Dstar);

% 4.1 Assembling and displaying matrix
abcd = [a,b;c,d];
disp('============================== Laminate Compliance Matrix (abcd) =============================');
disp(abcd);
disp('==============================================================================================');
disp(' ');

% 5.0 Effective Laminate Engineering Properties
% Average laminate stressess
Ex = 1/(Tt*a(1,1));
Ey = 1/(Tt*a(2,2));
Gxy = 1/(Tt*a(3,3));
vxy = -a(2,1)/a(1,1);
vyx = -a(1,2)/a(2,2);
n31 = a(1,3)/a(3,3);
n23 = a(3,2)/a(2,2);
% if sym == 0 %(is symmetrical)
%     n31 = 0;
%     n23 = 0;
% end

% 5.1 Displaying Engineering Properties
stiffness = [Ex,Ey,Gxy,vxy,vyx,n31,n23];
disp('=================== Average Effective Laminate Engineering Properties (Mpa) ==================');
disp(' ');
disp('         Ex           Ey           Gxy         vxy          vyx          n31          n23     ');
disp(stiffness);
disp('        [Pa]         [Pa]          [Pa] ');
disp('==============================================================================================');


% 6.0 Forces, Moments, and Stress Prediction
% Input Box
str = sprintf('Stiffness (ABD) and compliance (abcd) matrices are shown in the command \nwindow, along with the calculated effective laiminate engineering properties.');
str2 = sprintf('\n                                                                                                           Input forces and moments below. (Tension positive) ');
str = strcat(str, str2);

Nx = sprintf('\nNx (N/m) =');
str = strcat(str, Nx);
prompt = {str,'Ny =', 'Nxy =', 'Mx (Nm/m) =', 'My =', 'Mxy ='};
defaultanswer = {'100000','0','0','0','0','0'};
name = 'CLT - Force Input';
numlines = 1;
options.Resize = 'on';
options.WindowStyle = 'normal';
options.Interpreter = 'tex';
answer = inputdlg(prompt,name,numlines,defaultanswer,options);

% 6.0.1 Extracting data
for i = 1:1:6
F(i) = str2double(getfield(answer,{i}));
end

% 6.0.2 Force and moment input display table
disp(' ');
disp(' ');
disp('User input 2:');
disp(' ');
disp('---------------------------------------------------------------------------');
disp('                        Force and Moment Config                            ');
disp('---------------------------------------------------------------------------');
disp('       Nx(N/m)        Ny          Nxy      Mx(Nm/m)       My         Mxy   ');
disp(F);
disp('---------------------------------------------------------------------------');
disp('--> Forces loaded');
disp('--> Calculating strains and curvatures');

% 6.0.3 Calculating Strain and Curvature
ek = abcd*F.'; 
disp(' ');
disp(' ');
disp('====================== Strain and Curvature - Reference Surface ==============================');

disp('     Strain X     Strain Y    Strain XY  Curvature X  Curvature Y Curvature XY                ');
disp(ek');
disp('                                                                                       [mm/mm]');
disp('==============================================================================================');
disp(' ');


% 6.1 Stress Prediction for Each Layer
% 6.1.1 Global coordinate axes (x-y)
disp('--> Calculating stresses');
% Inverting h vector (calculate plys top to bottom)
h = -1.*h;
stressxy = zeros(lays*3,2);
p = 0;
for z = 1:3:lays*3
    p = p + 1;
    stressxy(z,1) = p;
    stressxy(z:z+2,2) = Qbar(z:z+2,1:3)*ek(1:3,1);
    % Surface stresses (global)
    stressxyt(z:z+2,2) = Qbar(z:z+2,1:3)*ek(1:3,1) + h(p)*Qbar(z:z+2,1:3)*ek(4:6,1); % Stress on top
    stressxyb(z:z+2,2) = Qbar(z:z+2,1:3)*ek(1:3,1) + h(p+1)*Qbar(z:z+2,1:3)*ek(4:6,1); % Stress on bottom
end

% 6.1.2 Local coordinate axes (x-y)
stress12 = zeros(lays*3,2);
for z = 1:3:lays*3
    stress12(z,1) = (z+2)/3;
    stress12(z:z+2,2) = Tmatrix(z:z+2,1:3)*stressxy(z:z+2,2);
    % Surface stresses (local)
    stress12t(z:z+2,2) = Tmatrix(z:z+2,1:3)*stressxyt(z:z+2,2);
    stress12b(z:z+2,2) = Tmatrix(z:z+2,1:3)*stressxyb(z:z+2,2);
end

% 6.2 Displaying Stresses (XY & 12)
plyt = zeros(lays,1);
for i = 1:1:lays
    plyt(i) = plytheta(i*3-2,1);
end
% 6.2.1 Table for GLOBAL Stress
gtable = zeros(lays,6);
for k = 1:1:lays
    gtable(k,1) = plyt(k);            % Ply
    gtable(k,2) = stressxyt(3*k-2,2); % Xt
    gtable(k,3) = stressxyb(3*k-2,2); % Xb
    gtable(k,4) = stressxyt(3*k-1,2); % Yt
    gtable(k,5) = stressxyb(3*k-1,2); % Yb
    gtable(k,6) = stressxyt(3*k,2);   % XYt
    gtable(k,7) = stressxyb(3*k,2);   % XYb
end
disp(' ');
disp('======================= Global (X-Y) Ply Stresses: Top and Bottom ============================');
disp('          Ply           Xt           Xb           Yt           Yb          XYt          XYb');
disp(gtable);
disp('                                                                                          [Pa]');
disp('==============================================================================================');

% 6.2.2 Table for LOCAL Stress
ltable = zeros(lays,6);
for k = 1:1:lays
    ltable(k,1) = plyt(k);            % Ply
    ltable(k,2) = stress12t(3*k-2,2); % Xt
    ltable(k,3) = stress12b(3*k-2,2); % Xb
    ltable(k,4) = stress12t(3*k-1,2); % Yt
    ltable(k,5) = stress12b(3*k-1,2); % Yb
    ltable(k,6) = stress12t(3*k,2);   % XYt
    ltable(k,7) = stress12b(3*k,2);   % XYb
end
disp(' ');
disp('======================= Local (1-2) Ply Stresses: Top and Bottom =============================');
disp('          Ply           1t           1b           2t           2b          12t          12b   ');
disp(ltable);
disp('                                                                                          [Pa]');
disp('==============================================================================================');




% 6.3 Configure S1, S2, S12 and ply reference matrices
for i = 1:3:lays*3
    c = (i+2)/3;
    S1(c,1) = stress12(i,2);
    S2(c,1) = stress12(i+1,2);
    S12(c,1) = stress12(i+2,2);
    plyref(c,1) = plymat(i,1);
end
stress = [S1,S2,S12];



% 6.4 Ultimate tensile strength matrix w.r.t. ply material (from Table 2.)
ults = zeros(lays,5);
for i = 1:1:lays
    ults(i,1:5) = prop2(plyref(i),2:6);
end


% 7.0 Check for Failure
% 7.1 Maximum Stress Theory
disp(' ');
disp(' ');
disp('---> Laminate Failure Predictions')
disp(' ');
disp('PART 1: Non-Interactive : Maximum Stress Theory');
disp(' ');
failtable = zeros(lays,11); 
flag = 0;
for i = 1:1:lays
    failtable(i,1) = i;
    failtable(i,2) = (plyt(i)); 
    if ltable(i,2) > 0 % 1T
        failtable(i,4) = 100*abs(ltable(i,2)/ults(i,1)); % percentage Ults tension
    else
        failtable(i,4) = 100*abs(ltable(i,2)/ults(i,3)); % percentage Ults compression
    end
    if ltable(i,3) > 0 % 1B
        failtable(i,6) = 100*abs(ltable(i,3)/ults(i,1)); % percentage Ults tension
    else
        failtable(i,6) = 100*abs(ltable(i,3)/ults(i,3)); % percentage Ults compression
    end
    if ltable(i,4) > 0 % 2T
        failtable(i,8) = 100*abs(ltable(i,4)/ults(i,2)); % percentage Ults tension
    else
        failtable(i,8) = 100*abs(ltable(i,4)/ults(i,4)); % percentage Ults compression
    end
    if ltable(i,4) > 0 % 2B
        failtable(i,10) = 100*abs(ltable(i,5)/ults(i,2)); % percentage Ults tension
    else
        failtable(i,10) = 100*abs(ltable(i,5)/ults(i,4)); % percentage Ults compression
    end
    failtable(i,12) = 100*abs(ltable(i,6)/ults(i,5)); % percentage Ults SHEAR
    
    % Key:
    % 0 = OK in Tension, 1 = OK in Compression, 2 = OK in Shear
    % 7 = FAILURE in Tension, 8 = FAILURE in Comp, 9 = FAILURE in Shear
    
    c = 0;
    for w = 3:2:8
        c = c + 1;
        if ltable(i,c+1) > 0 % IN TENSION
            if failtable(i,w+1) > 100 % Sigma 1t - Tension
                failtable(i,w) = 7; %('FAILURE! (T))');
                flag = 1;
            else 
                failtable(i,w) = 0; %('Ok in T');
            end
        else % IN COMPRESSION
            if failtable(i,w+1) > 100 % Sigma 1t - Compression
                failtable(i,w) = 8; %('FAILURE! (C)');
                flag = 1;
            else 
                failtable(i,w) = 1; %('Ok in C');
            end
        end   
    end
    if failtable(i,12) > 100 % Tau12
        failtable(i,11) = 9; % ('Failure in Shear');
        flag = 1;
    else
        failtable(i,11) = 2; % ('Ok in Shear');
    end
end

% 7.1.1 Display Failtable

disp('========================= Non-Interactive Failure Prediction : Maximum Stress Theory =========================================================================');
disp('          Ply          Dir           1t         %Ult           1b         %Ult           2t         %Ult           2b         %Ult           12         %Ult');
disp(failtable);
disp('--------------------------------------------------------------------------------------------------------------------------------------------------------------')
disp('          Key: 0 = OK in Tension, 1 = OK in Compression, 2 = OK in Shear');
disp('               7 = FAILURE in Tension, 8 = FAILURE in Compression, 9 = FAILURE in Shear');
disp('==============================================================================================================================================================');

% 7.1.2 Identifying  max %Ult in each field - then displaying the
% corresponding (max) ply number
[v(1,1),v(1,2)] = max(failtable(:,3));
[v(2,1),v(2,2)] = max(failtable(:,5));
[v(3,1),v(3,2)] = max(failtable(:,7));
[v(4,1),v(4,2)] = max(failtable(:,9));
[v(5,1),v(5,2)] = max(failtable(:,11));
[vf,p]  = max(v(:,1));
pf = v(p,2);

if flag == 0
    str = sprintf('--> Material OK under loading! First Ply Failure (FPF) will occur in ply %g.', pf);
    disp(str)
else
    str = sprintf('--> Material FAILURE under loading! First Ply Failure (FPF) occurs in ply %g. (See table above)', pf);
    disp(str);
end



% 7.2 Tsai-Hill Theory
% Page 87 C. MACRO
disp(' ');
disp(' ');
disp('PART 2. Interactive : Tsai-Hill Theory');
flag = 0;
failtableTH = zeros(lays,4); 

for i = 1:1:lays
    failtableTH(i,1) = i;
    failtableTH(i,2) = plyt(i);
    
    % Tsai-Hill for TOP SURFACE
    if ltable(i,2) > 0
        s1term = ltable(i,2)/ults(i,1); % Tension
    else
        s1term = -ltable(i,2)/ults(i,3); % Compression
    end
    if ltable(i,4) > 0
        s2term = ltable(i,4)/ults(i,2); % Tension 
    else
        s2term = -ltable(i,4)/ults(i,4); % Compression
    end
    t12term = abs(ltable(i,6)/ults(i,5));
    
    tsai = s1term^2 + s2term^2 + t12term^2 - s1term*s2term;
    if tsai > 1
        flag = 1;
        failtableTH(i,3) = 9; % Failure
    else
        failtableTH(i,3) = 0;
    end
    
    % Tsai-Hill for BOTTOM SURFACE
    if ltable(i,3) > 0
        s1term = ltable(i,3)/ults(i,1); % Tension
    else
        s1term = -ltable(i,3)/ults(i,3); % Compression
    end
    if ltable(i,4) > 0
        s2term = ltable(i,5)/ults(i,2); % Tension 
    else
        s2term = -ltable(i,5)/ults(i,4); % Compression
    end
    t12term = abs(ltable(i,7)/ults(i,5));
    
    tsai = s1term^2 + s2term^2 + t12term^2 - s1term*s2term;
    if tsai > 1
        flag = 1;
        failtableTH(i,4) = 9; % Failure
    else
        failtableTH(i,4) = 0; % OK!
    end
end

% Display Failtable - Tsai-Hill

disp(' ');
disp('============ Interactive Failure Prediction : Tsai-Hill Theory ===============');
disp('                 Surface   ');
disp('   Ply   Dir   Top  Bottom ');
disp(failtableTH);
disp('------------------------------------------------------------------------------');
disp('   Key: 0 = Material OK, 9 = Material FAILURE');
disp('  (Note: the Tsai-Hill criterion only indicates whether or not');
disp('   a ply will fail, it does not indicate the mode of failure.)');
disp('==============================================================================');


if flag == 0
    disp('--> Material OK under loading! (Tsai-Hill criterion)');
else
    disp('--> Material FAILURE under loading! (Tsai-Hill criterion) See table above.');
end

% End program






