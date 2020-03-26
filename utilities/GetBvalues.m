function [bval]= GetBvalues(protocol)
%
% camino.m--------------------------------------------------------------
% Computes b-values
%
% [bval]= GetBvalues(protocol)
% 
% Description: Returns an array of b values, one for each measurement defined
% in the list of measurements in the protocol.
%
% Paramaters: 
% b - list of b-values
% protocol - diffusion acquisition protocol 
%
%------------------------------------------------------------------------
% This file is part of the camino.m toolbox.
% Copyright (c) 2015, UCL Microstructure Imaging Group (MIG), All rights reserved.
% Distributed under the Modified BSD Licence (see: LICENSE.pdf).
%
% Author:
%   Daniel C Alexander (d.alexander@ucl.ac.uk)
%   Andrada Ianus (a.ianus.11@ucl.ac.uk)
%



GAMMA = 2.675987E8;
if(strcmp(protocol.pulseseq, 'PGSE') || strcmp(protocol.pulseseq,'STEAM'))
    modQ = GAMMA*protocol.smalldel.*protocol.G;
    diffTime = protocol.delta - protocol.smalldel/3;
    bval = diffTime.*modQ.^2;
elseif strcmp(protocol.pulseseq,'DSE')
    G = protocol.G; 
    t1 = protocol.t1; t2 = protocol.t2; t3 = protocol.t3;
    delta1 = protocol.delta1; delta2 = protocol.delta2; delta3 = protocol.delta3;
    
    bval = (GAMMA*G).^2;
    bval = bval.*(-t1.*(delta1.^2) + t3.*(delta1.^2) - (delta1.^3)./3 ...
                - 2*t2.*delta1.*delta2 + 2*t3.*delta1.*delta2 ...
                + (delta1.^2).*delta2 - t2.*(delta2.^2) + t3.*(delta2.^2) ...
                - (delta2.^3)./3 + 2*t2.*delta1.*delta3 - 2*t3.*delta1.*delta3 ...
                - (delta1.^2).*delta3 + 2*t2.*delta2.*delta3 ...
                - 2*t3.*delta2.*delta3 + (delta2.^2).*delta3 ...
                - t2.*(delta3.^2) + t3.*(delta3.^2) + 2*delta1.*(delta3.^2)  ...
                + delta2.*(delta3.^2) - (delta3.^3));

elseif(strcmp(protocol.pulseseq, 'OGSE'))
    if ~isfield(protocol,'mirror') || protocol.mirror == 0
         G = protocol.G;
         omega = protocol.omega;
         delta = protocol.delta;
         smalldel = protocol.smalldel;
        if ~isfield(protocol,'phase')           
            bval = (GAMMA*G).^2.*(delta.*omega.*(3-4*cos(omega.*smalldel)+cos(2*omega.*smalldel)) +...
            omega.*smalldel.*(2 + 4*cos(omega.*smalldel)) - 4*sin(omega.*smalldel) - sin(2*omega.*smalldel))./(2*omega.^3);
        else
            phase = protocol.phase;
           bval = (GAMMA*G).^2./(2*omega.^3).*(2.*(delta-smalldel).*omega.*(cos(phase)-cos(smalldel.*omega-phase)).^2 +...
        4.*smalldel.*omega - 4.*sin(smalldel.*omega)- sin(2.*phase)+smalldel.*omega.*(cos(2*phase)+...
            cos(2.*smalldel.*omega-2.*phase))-sin(2.*smalldel.*omega-2.*phase));
        end
    else
        error('Not yet implemented')
    end
      

elseif(strcmp(protocol.pulseseq, 'SWOGSE'))     % not reflected
    G = protocol.G;
    omega = protocol.omega;    
    niu = omega/(2*pi());
    smalldel = protocol.smalldel;
    delta = protocol.delta;
    
    if ~isfield(protocol,'mirror') || protocol.mirror == 0
        if ~isfield(protocol,'phase') 
            NT = floor(2.*smalldel.*niu+0.00000000001);
            bval = GAMMA.^2.*G.^2.*(1./(48.*niu.^3)).*((2.*NT.^3 +... 
             3.*NT.^2.*(1 + (-1).^NT - 4.*niu.*smalldel) - 4.*niu.^2.*smalldel.^2.*...
             (-3 + 3.*(-1).^NT + 4.*niu.*smalldel) +3.*delta.*niu.*(-1 + (-1).^NT - ...
             2.*NT + 4.*niu.*smalldel).^2 + NT.*(1 + 3.*(-1).^NT - 12.*niu.*smalldel + ...
               24.*niu.^2.*smalldel.^2)));
        else
             phase = protocol.phase;
             for i = 1:length(omega)
                if omega(i)<pi/smalldel(i);
                    omega(i) = pi/smalldel(i);
                    phase(i) = 0;
                end
             end
            phase = mod(phase,2*pi);
            phase(phase>pi) = phase(phase>pi)-2*pi;
            phase(phase<0) = pi-abs(phase(phase<0));
            phdelay = phase ./(2 *pi()* niu);

            NT = floor(2.*(smalldel-phdelay).*niu+0.00000000001);
            sgn = (-1).^NT;

            bval = GAMMA.^2.*G.^2.*(1./3.*(smalldel-NT./2./niu-phdelay).^3 + (delta-smalldel).*...
            (sgn.*(smalldel- (0.5.*(1-sgn)+NT)./2./niu-phdelay)-phdelay).^2  +(phdelay.^3)./3 +...
            sgn.*((-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay)).^3)./192./(niu.^3) +...   
             NT./96./(niu.^3).*(8+12.*NT.*(1+NT)- 24.*smalldel.*niu.*(1 + 2.*NT) +48.*smalldel.^2.*niu.^2+...
            48.*NT.*niu.*phdelay - 96.*niu.^2.*phdelay.*(smalldel-phdelay)) +...
            1./3.*(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3 +...     
            1/3.*sgn.*((-1+sgn+4.*niu.*phdelay).^3./64./niu.^3-(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3));

        end
    else
         if ~isfield(protocol,'phase') 
             NT = floor(2.*smalldel.*niu+0.00000000001);

            bval = GAMMA.^2.*G.^2.*((delta-smalldel).*(smalldel-(0.5.*(1-(-1).^NT)+NT)./(2.*niu)).^2 + NT./(12.*niu.^3) +...
            1./(192.*niu.^3).*(-(-1+(-1).^NT).^3+(-1+(-1).^NT-2.*NT+4.*smalldel.*niu).^3) - ...
            1./(96.*niu.^3).*(NT-2.*smalldel.*niu).*(3.*(-1+(-1).^NT).^2+4.*NT.^2+12.*smalldel.*niu.*(-1+(-1).^NT)+...
            16.*smalldel.^2.*niu.^2-2.*NT.*(-3+3.*(-1).^NT+8.*smalldel.*niu)));
         else
              phase = protocol.phase;
             for i = 1:length(omega)
                if omega(i)<pi/smalldel(i);
                    omega(i) = pi/smalldel(i);
                    phase(i) = 0;
                end
             end
            phase = mod(phase,2*pi);
            phase(phase>pi) = phase(phase>pi)-2*pi;
            phase(phase<0) = pi-abs(phase(phase<0));
            phdelay = phase ./(2 *pi()* niu);

            NT = floor(2.*(smalldel-phdelay).*niu+0.00000000001);
            sgn = (-1).^NT;

            bval = G.^2.*GAMMA.^2.*((delta-smalldel).*(sgn.*(smalldel-(0.5.*(1-sgn)+NT)./(2.*niu)-phdelay)-phdelay).^2 +...
            2.*phdelay.^3./3+ NT.*(1-6.*niu.*phdelay+12.*niu.^2.*phdelay.^2)./12./niu.^3+...
            sgn./3.*((phdelay-sgn./4./niu.*(sgn-1)).^3-2.*(phdelay-sgn./4./niu.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))).^3)+...
            sgn./3.*((sgn-1+4.*niu.*phdelay).^3./64./niu.^3));
         end
    end
   

elseif (strcmp(protocol.pulseseq, 'TWOGSE'))
    delta = protocol.delta;
    smalldel = protocol.smalldel;
    G = protocol.G;
    omega = protocol.omega;    
    niu = omega/(2*pi());
    slew_rate = protocol.slew_rate;
    NT = floor(2.*smalldel.*niu+0.00000000001); 
    rt = G./slew_rate;
%     if max(protocol.smalldel - NT./2./niu) > 1E-5 % calculate bvalue from waveform
%         tau = 1E-4;
%         protocol.tau = tau;
%         wf = wave_form(protocol);
%         Fx = cumsum(wf(:,1:3:end).*tau,2);
%         Fy = cumsum(wf(:,2:3:end).*tau,2);
%         Fz = cumsum(wf(:,3:3:end).*tau,2);
%         bval=sum((Fx.^2+Fy.^2+Fz.^2)*tau,2);   
%         
%     else        
        bval = GAMMA.^2.*G.^2.*((delta-smalldel)./4.*(1-(-1).^NT).^2.*(1./2./niu-rt).^2+...
        smalldel./(240.*niu.^2).*(40-120.*rt.*niu-40.*rt.^2.*niu.^2+256.*rt.^3.*niu.^3));
%     end


elseif strcmp(protocol.pulseseq, 'SWOGSE_3D') 
    tau=protocol.tau;         
    K = floor((protocol.smalldel+1E-10)./tau)+1;
    dstmp = floor((protocol.delta-1E-10)./protocol.tau)-floor((protocol.smalldel+1E-10)./protocol.tau);
    if ~isfield(protocol,'angle') || protocol.angle == 4
          Gx = protocol.Gx';
          Gy = protocol.Gy';
          Gz = protocol.Gz';
        M = size(Gx,1);
        Bval=zeros(1,M);      
         if ~isfield(protocol,'phix') && ~isfield(protocol,'phiy') && ~isfield(protocol,'phiz')
            phix = zeros(1,M); phiy = zeros(1,M); phiz = zeros(1,M);
        else
            phix = protocol.phix; phiy = protocol.phiy;  phiz = protocol.phiz;
        end

            for m=1:M                 

                    time_vec = tau*(0:K(m)-1);
                    sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-phix(m))./pi-1E-10);
                    sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-phiy(m))./pi-1E-10);
                    sign_vec3 = (-1).^floor((protocol.omegaz(m)*time_vec-phiz(m))./pi-1E-10);
                    vec1 = zeros(size(sign_vec1)); vec2 = zeros(size(sign_vec2));  vec3 = zeros(size(sign_vec3)); 
                    vec1(sign_vec1>0) = Gx(m); vec1(sign_vec1<0) = -Gx(m); vec1(1)=0; vec1(end)=0;
                    vec2(sign_vec2>0) = Gy(m); vec2(sign_vec2<0) = -Gy(m); vec2(1)=0; vec2(end)=0;
                    vec3(sign_vec3>0) = Gz(m); vec3(sign_vec3<0) = -Gz(m); vec3(1)=0; vec3(end)=0;


                      if ~isfield(protocol,'mirror') || protocol.mirror == 0
                          Gx_vec = [vec1 zeros(1,dstmp(m)) -(vec1)]; 
                        Gy_vec = [vec2 zeros(1,dstmp(m)) -(vec2)]; 
                        Gz_vec = [vec3 zeros(1,dstmp(m)) -(vec3)]; 
                      else
                        Gx_vec = [vec1 zeros(1,dstmp(m)) -fliplr(vec1)]; 
                        Gy_vec = [vec2 zeros(1,dstmp(m)) -fliplr(vec2)]; 
                        Gz_vec = [vec3 zeros(1,dstmp(m)) -fliplr(vec3)];                      
                        
                      end

                    Fx = cumsum(Gx_vec)*tau;
                    Fy = cumsum(Gy_vec)*tau;
                    Fz = cumsum(Gz_vec)*tau;
                    Bval(m)=sum((Fx.^2+Fy.^2+Fz.^2)*tau);              



            end
            bval=GAMMA^2*Bval;
         
    
    elseif protocol.angle == 1 % gradient only along x
        Gx = protocol.Gx';
        M = size(Gx,1);   
        Bval=zeros(1,M);    
        if ~isfield(protocol,'phix')
            phix = zeros(1,M);
        else
            phix = protocol.phix;
        end
      

        for m=1:M
                time_vec = tau*(0:K(m)-1);
                sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-phix(m))./pi-1E-10);
                vec1 = zeros(size(sign_vec1));
                vec1(sign_vec1>0) = Gx(m); vec1(sign_vec1<0) = -Gx(m); vec1(1)=0; vec1(end)=0;
                 if ~isfield(protocol,'mirror') || protocol.mirror == 0   
                     Gx_vec = [vec1 zeros(1,dstmp(m)) -(vec1)];
                 else                
                     Gx_vec = [vec1 zeros(1,dstmp(m)) -fliplr(vec1)];                   
                 end

                Fx = cumsum(Gx_vec)*tau;
                Bval(m)=sum((Fx.^2)*tau);        

        end
    bval=GAMMA^2*Bval;     
    elseif protocol.angle == 2
      Gx = protocol.Gx';
      Gy = protocol.Gy';
        M = size(Gx,1);
        Bval=zeros(1,M);    
        if ~isfield(protocol,'phix') && ~isfield(protocol,'phiy')
            phix = zeros(1,M); phiy = zeros(1,M);
        else
            phix = protocol.phix; phiy = protocol.phiy; 
        end

            for m=1:M                 

                    time_vec = tau*(0:K(m)-1);
                    sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-phix(m))./pi-1E-10);
                    sign_vec2 = (-1).^floor((protocol.omegay(m)*time_vec-phiy(m))./pi-1E-10);
                    vec1 = zeros(size(sign_vec1)); vec2 = zeros(size(sign_vec2)); 
                    vec1(sign_vec1>0) = Gx(m); vec1(sign_vec1<0) = -Gx(m); vec1(1)=0; vec1(end)=0;
                    vec2(sign_vec2>0) = Gy(m); vec2(sign_vec2<0) = -Gy(m); vec2(1)=0; vec2(end)=0;


                      if ~isfield(protocol,'mirror') || protocol.mirror == 0    
                        Gx_vec = [vec1 zeros(1,dstmp(m)) -(vec1)]; 
                        Gy_vec = [vec2 zeros(1,dstmp(m)) -(vec2)]; 
                      else
                         Gx_vec = [vec1 zeros(1,dstmp(m)) -fliplr(vec1)]; 
                        Gy_vec = [vec2 zeros(1,dstmp(m)) -fliplr(vec2)]; 
                      end

                    Fx = cumsum(Gx_vec)*tau;
                    Fy = cumsum(Gy_vec)*tau;
                    Bval(m)=sum((Fx.^2+Fy.^2)*tau);                 



            end
            bval=GAMMA^2*Bval;
            
    elseif protocol.angle == 3
          Gx = protocol.Gx';
          Gz = protocol.Gz';
        M = size(Gx,1);
        Bval=zeros(1,M);  
         if ~isfield(protocol,'phix') && ~isfield(protocol,'phiz')
            phix = zeros(1,M); phiz = zeros(1,M);
        else
            phix = protocol.phix; phiz = protocol.phiz; 
        end

            for m=1:M                 

                    time_vec = tau*(0:K(m)-1);
                    sign_vec1 = (-1).^floor((protocol.omegax(m)*time_vec-phix(m))./pi-1E-10);
                    sign_vec2 = (-1).^floor((protocol.omegaz(m)*time_vec-phiz(m))./pi-1E-10);
                    vec1 = zeros(size(sign_vec1)); vec2 = zeros(size(sign_vec2)); 
                    vec1(sign_vec1>0) = Gx(m); vec1(sign_vec1<0) = -Gx(m); vec1(1)=0; vec1(end) = 0;
                    vec2(sign_vec2>0) = Gz(m); vec2(sign_vec2<0) = -Gz(m); vec2(1)=0; vec2(end) = 0;


                      if  ~isfield(protocol,'mirror') || protocol.mirror == 0
                     
                        Gx_vec = [vec1 zeros(1,dstmp(m)) -(vec1)]; 
                        Gz_vec = [vec2 zeros(1,dstmp(m)) -(vec2)];
                      else
                        Gx_vec = [vec1 zeros(1,dstmp(m)) -fliplr(vec1)]; 
                        Gz_vec = [vec2 zeros(1,dstmp(m)) -fliplr(vec2)]; 
                      
                      end

                    Fx = cumsum(Gx_vec)*tau;
                    Fz = cumsum(Gz_vec)*tau;
                    Bval(m)=sum((Fx.^2+Fz.^2)*tau);                


            end
            bval=GAMMA^2*Bval;
           
    end
elseif strcmp(protocol.pulseseq,'dPGSE')
    if ~isfield(protocol,'slew_rate')
    
         if sum(protocol.tm >= protocol.smalldel) == length(protocol.tm)
            b1 = GAMMA.^2.*protocol.G1.^2.*protocol.smalldel.^2.*(protocol.delta-protocol.smalldel/3);   
            b2 = GAMMA.^2.*protocol.G2.^2.*protocol.smalldel.^2.*(protocol.delta-protocol.smalldel/3);
            bval = b1 + b2;
         else
             error('bvalue not defined yet for tm < smalldel');
         end
    else
        slew_rate = protocol.slew_rate;
        rt1 = protocol.G1./slew_rate;
        rt2 = protocol.G2./slew_rate;  
        delta1 = protocol.delta;
        smalldel1 = protocol.smalldel;
        delta2 = protocol.delta;
        smalldel2 = protocol.smalldel;        
        if sum(protocol.tm >= protocol.smalldel) == length(protocol.tm)
        bval = 2*GAMMA.^2.*(protocol.G1.^2.*smalldel1.^3./15.*(5-15.*rt1./smalldel1./2 - ...
            5.*rt1.^2./smalldel1./4 + 4.*rt1.^3./smalldel1.^3)) +...  
            GAMMA.^2.*protocol.G1.^2.*(delta1 - smalldel1).*(smalldel1-rt1).^2 + ...
            2*GAMMA.^2.*(protocol.G2.^2.*smalldel2.^3./15.*(5-15.*rt2./smalldel2./2 - ...
            5.*rt2.^2./smalldel2./4 + 4.*rt2.^3./smalldel2.^3)) +...  
            GAMMA.^2.*protocol.G2.^2.*(delta2 - smalldel2).*(smalldel2-rt2).^2;
         else
             error('bavlue not defined yet for tm < smalldel');
         end
    end
   
elseif strcmp(protocol.pulseseq,'dSWOGSE')
    G1 = protocol.G1;
    G2 = protocol.G2;
    omega = protocol.omega;    
    niu = omega/(2*pi());
    smalldel = protocol.smalldel;
    delta = protocol.delta;
    
    if ~isfield(protocol,'mirror') || protocol.mirror == 0
        if ~isfield(protocol,'phase') 
            NT = floor(2.*smalldel.*niu+0.00000000001);
            b1 = GAMMA.^2.*G1.^2.*(1./(48.*niu.^3)).*((2.*NT.^3 +... 
             3.*NT.^2.*(1 + (-1).^NT - 4.*niu.*smalldel) - 4.*niu.^2.*smalldel.^2.*...
             (-3 + 3.*(-1).^NT + 4.*niu.*smalldel) +3.*delta.*niu.*(-1 + (-1).^NT - ...
             2.*NT + 4.*niu.*smalldel).^2 + NT.*(1 + 3.*(-1).^NT - 12.*niu.*smalldel + ...
               24.*niu.^2.*smalldel.^2)));
            b2 = GAMMA.^2.*G2.^2.*(1./(48.*niu.^3)).*((2.*NT.^3 +... 
             3.*NT.^2.*(1 + (-1).^NT - 4.*niu.*smalldel) - 4.*niu.^2.*smalldel.^2.*...
             (-3 + 3.*(-1).^NT + 4.*niu.*smalldel) +3.*delta.*niu.*(-1 + (-1).^NT - ...
             2.*NT + 4.*niu.*smalldel).^2 + NT.*(1 + 3.*(-1).^NT - 12.*niu.*smalldel + ...
               24.*niu.^2.*smalldel.^2)));
        else
             phase = protocol.phase;
             for i = 1:length(omega)
                if omega(i)<pi/smalldel(i);
                    omega(i) = pi/smalldel(i);
                    phase(i) = 0;
                end
             end
            phase = mod(phase,2*pi);
            phase(phase>pi) = phase(phase>pi)-2*pi;
            phase(phase<0) = pi-abs(phase(phase<0));
            phdelay = phase ./(2 *pi()* niu);

            NT = floor(2.*(smalldel-phdelay).*niu+0.00000000001);
            sgn = (-1).^NT;

            b1 = GAMMA.^2.*G1.^2.*(1./3.*(smalldel-NT./2./niu-phdelay).^3 + (delta-smalldel).*...
                (sgn.*(smalldel- (0.5.*(1-sgn)+NT)./2./niu-phdelay)-phdelay).^2  +(phdelay.^3)./3 +...
                sgn.*((-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay)).^3)./192./(niu.^3) +...   
                 NT./96./(niu.^3).*(8+12.*NT.*(1+NT)- 24.*smalldel.*niu.*(1 + 2.*NT) +48.*smalldel.^2.*niu.^2+...
                48.*NT.*niu.*phdelay - 96.*niu.^2.*phdelay.*(smalldel-phdelay)) +...
                1./3.*(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3 +...     
                1/3.*sgn.*((-1+sgn+4.*niu.*phdelay).^3./64./niu.^3-(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3));
            b2 = GAMMA.^2.*G2.^2.*(1./3.*(smalldel-NT./2./niu-phdelay).^3 + (delta-smalldel).*...
                (sgn.*(smalldel- (0.5.*(1-sgn)+NT)./2./niu-phdelay)-phdelay).^2  +(phdelay.^3)./3 +...
                sgn.*((-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay)).^3)./192./(niu.^3) +...   
                 NT./96./(niu.^3).*(8+12.*NT.*(1+NT)- 24.*smalldel.*niu.*(1 + 2.*NT) +48.*smalldel.^2.*niu.^2+...
                48.*NT.*niu.*phdelay - 96.*niu.^2.*phdelay.*(smalldel-phdelay)) +...
                1./3.*(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3 +...     
                1/3.*sgn.*((-1+sgn+4.*niu.*phdelay).^3./64./niu.^3-(phdelay-sgn.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))./4./niu).^3));

        end
    else
         if ~isfield(protocol,'phase') 
             NT = floor(2.*smalldel.*niu+0.00000000001);

            b1 = GAMMA.^2.*G1.^2.*((delta-smalldel).*(smalldel-(0.5.*(1-(-1).^NT)+NT)./(2.*niu)).^2 + NT./(12.*niu.^3) +...
                1./(192.*niu.^3).*(-(-1+(-1).^NT).^3+(-1+(-1).^NT-2.*NT+4.*smalldel.*niu).^3) - ...
                1./(96.*niu.^3).*(NT-2.*smalldel.*niu).*(3.*(-1+(-1).^NT).^2+4.*NT.^2+12.*smalldel.*niu.*(-1+(-1).^NT)+...
                16.*smalldel.^2.*niu.^2-2.*NT.*(-3+3.*(-1).^NT+8.*smalldel.*niu)));
            b2 = GAMMA.^2.*G2.^2.*((delta-smalldel).*(smalldel-(0.5.*(1-(-1).^NT)+NT)./(2.*niu)).^2 + NT./(12.*niu.^3) +...
                1./(192.*niu.^3).*(-(-1+(-1).^NT).^3+(-1+(-1).^NT-2.*NT+4.*smalldel.*niu).^3) - ...
                1./(96.*niu.^3).*(NT-2.*smalldel.*niu).*(3.*(-1+(-1).^NT).^2+4.*NT.^2+12.*smalldel.*niu.*(-1+(-1).^NT)+...
                16.*smalldel.^2.*niu.^2-2.*NT.*(-3+3.*(-1).^NT+8.*smalldel.*niu)));
         else
              phase = protocol.phase;
             for i = 1:length(omega)
                if omega(i)<pi/smalldel(i);
                    omega(i) = pi/smalldel(i);
                    phase(i) = 0;
                end
             end
            phase = mod(phase,2*pi);
            phase(phase>pi) = phase(phase>pi)-2*pi;
            phase(phase<0) = pi-abs(phase(phase<0));
            phdelay = phase ./(2 *pi()* niu);

            NT = floor(2.*(smalldel-phdelay).*niu+0.00000000001);
            sgn = (-1).^NT;

            b1 = G1.^2.*GAMMA.^2.*((delta-smalldel).*(sgn.*(smalldel-(0.5.*(1-sgn)+NT)./(2.*niu)-phdelay)-phdelay).^2 +...
                2.*phdelay.^3./3+ NT.*(1-6.*niu.*phdelay+12.*niu.^2.*phdelay.^2)./12./niu.^3+...
                sgn./3.*((phdelay-sgn./4./niu.*(sgn-1)).^3-2.*(phdelay-sgn./4./niu.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))).^3)+...
                sgn./3.*((sgn-1+4.*niu.*phdelay).^3./64./niu.^3));
            b2 = G2.^2.*GAMMA.^2.*((delta-smalldel).*(sgn.*(smalldel-(0.5.*(1-sgn)+NT)./(2.*niu)-phdelay)-phdelay).^2 +...
                2.*phdelay.^3./3+ NT.*(1-6.*niu.*phdelay+12.*niu.^2.*phdelay.^2)./12./niu.^3+...
                sgn./3.*((phdelay-sgn./4./niu.*(sgn-1)).^3-2.*(phdelay-sgn./4./niu.*(-1+sgn-2.*NT+4.*niu.*(smalldel-phdelay))).^3)+...
                sgn./3.*((sgn-1+4.*niu.*phdelay).^3./64./niu.^3));
         end
    end
    bval = b1+b2;
elseif strcmp(protocol.pulseseq,'isoPGSE') 
     bval = GAMMA.^2.*6*protocol.G.^2.*(protocol.smalldel./6).^2.*(2*protocol.smalldel/6/3);
elseif strcmp(protocol.pulseseq,'DODE')
    G1 = protocol.G1;
    G2 = protocol.G2;
    N = protocol.Nosc; % the number of periods
    smalldel = protocol.smalldel;

    tm = protocol.tm;
    
    if ~isfield(protocol,'slew_rate')
    
        if ~isfield(protocol,'phase') || max(protocol.phase) == 0      


        bval = GAMMA.^2.*(G1.^2+G2.^2)./2.*smalldel.^3./6./N.^2;
        elseif max(protocol.phase) == pi/2 || max(protocol.phase) == -pi/2
            bval = GAMMA.^2.*(G1.^2+G2.^2)./2.*smalldel.^3./24./N.^2;

        else error('DODE defined only for phase= 0 or pi/2')
        end
    else
     slew_rate = protocol.slew_rate;
    rt1 = G1./slew_rate;
    rt2 = G2./slew_rate;       
         if ~isfield(protocol,'phase') || max(protocol.phase) == 0     
               bval = GAMMA.^2.*(G1.^2.*smalldel.^3./15./N.^2./4.*(5-15.*rt1.*N./smalldel - ...
                    5.*rt1.^2.*N.^2./smalldel+4.*rt1.^3.*8.*N.^3./smalldel.^3) +...    
                     G2.^2.*smalldel.^3./15./N.^2./4.*(5-15.*rt2.*N./smalldel - ...
                    5.*rt2.^2.*N.^2./smalldel+4.*rt2.^3.*8.*N.^3./smalldel.^3));
             
         elseif max(protocol.phase) == pi/2 || max(protocol.phase) == -pi/2
         
               bval = GAMMA.^2.*(G1.^2.*smalldel.^3./12./4./N.^2.*((1+16.*2.*N).*rt1.^3.*4.*N.^2./5./smalldel.^3-...
                    4.*rt1.^2.*4.*N.^2./smalldel.^2+1) +...
                     G2.^2.*smalldel.^3./12./4./N.^2.*((1+16.*2.*N).*rt2.^3.*4.*N.^2./5./smalldel.^3-...
                    4.*rt2.^2.*4.*N.^2./smalldel.^2+1));
        else error('DODE defined only for phase= 0 or pi/2')
        end
    end
     
   
elseif strcmp(protocol.pulseseq,'GEN') 
     wf =protocol.G;
     tau = protocol.tau;
    Fx = cumsum(wf(:,1:3:end).*tau,2);
    Fy = cumsum(wf(:,2:3:end).*tau,2);
    Fz = cumsum(wf(:,3:3:end).*tau,2);
    bval=sum((Fx.^2+Fy.^2+Fz.^2)*tau,2);   
    bval = GAMMA^2.*bval';
else
error('the protocol does not contain enough information to generate the waveform');

end

end
