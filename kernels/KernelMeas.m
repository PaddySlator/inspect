function S = KernelMeas(kernel,gradechoinv)

%and returns a function handle for that kernel

%protocol gives the experimental  parameters  which  are  varied  to  yield  contrast in intrinsic MR properties w.



% g = protocol.graddirs;
% b = protocol.b;
% TE = protocol.TE;
% TI = protocol.TI;
% TR = protocol.TR;


fun = str2func(['Kernel' kernel.name]);

S = fun(kernel.params,gradechoinv);




% switch options.kernel
%     case 'ZeppelinT2T1'
%         fkernel = @(params,protocol) ...
%             Zeppelin(params,protocol) .* ...
%             t2decay(params,protocol) .* ...
%             t1inv(params,protocol);
%         
%     case 'BallT2T1'
%         %         fkernel = @(g,b,TE,TI,TR,d,T2,T1) ...
%         %             exp(-b*d) * ...
%         %             exp(-TE/T2) * ...
%         %             abs(1 - 2 * exp( -(TI + TE)/ T1) + exp(-TR/T1) );
%         
%         fkernel = @(params,protocol) ...                        
%             Ball(params,protocol) .* ...
%             t2decay(params,protocol) .* ...
%             t1inv(params,protocol);
%         
%     case 'BallBall'
%          fkernel = @(params,protocol) ...
%             Ball(params,protocol) .* ...
%             Ball(params,protocol);
%         
%     case 'BallT2'        
%         fkernel = @(params,protocol) ...
%             Ball(params,protocol) .* ...
%             t2decay(params,protocol);
%         
%     otherwise
%         error('Unknown kernel')
        

end
        








