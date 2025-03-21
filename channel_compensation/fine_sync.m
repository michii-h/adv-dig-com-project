function [Rlk_out] = fine_sync(Rlk_in, Plk, iNfft, SwitchDemoSync)
    % FINE_SYNC Linear phase correction using pilot tones
    % Performs fine synchronization by estimating and correcting linear phase errors
    %
    % Inputs:
    %   Rlk_in - Received OFDM frame in frequency domain
    %   Plk - Pilot pattern matrix (same size as Rlk_in)
    %   iNfft - FFT length
    %   SwitchDemoSync - Enable debug plots
    %
    % Outputs:
    %   Rlk_out - Phase-corrected OFDM frame

    % Copy input to output
    Rlk_out = Rlk_in;

    % Calculate phase differences between received signal and pilots
    mPhase = angle(conj(Rlk_in).*Plk);

    % Display 3D mesh of unwrapped phase if debug is enabled
    if SwitchDemoSync
        figure(301)
        mesh([-iNfft/2:iNfft/2-1],[1:size(Plk,1)],unwrap(mPhase))
        ylabel('Symbol l')
        xlabel('Subchannel k')
        zlabel('Phase \Phi(l,k)')
    end

    % Process each OFDM symbol individually
    for l = 1:size(Plk,1)
        % Display phase plot if debug is enabled
        if SwitchDemoSync && l == 1
            figure(302)
            stem([-iNfft/2:iNfft/2-1],unwrap(mPhase(l,:)))
            xlabel('Subchannel k')
            ylabel(['Phase \Phi(' num2str(l) ',k)'])
        end

        % Linear Regression of Phases
        kPilots = [find(Plk(l,:) ~= 0)-iNfft/2-1].';
        V = [kPilots ones(size(kPilots))];
        vPhi = unwrap(mPhase(l,Plk(l,:) ~= 0)).';
        m = V\vPhi;

        % Display linear fit if debug is enabled
        if SwitchDemoSync && l == 1
            hold on
            kAxis = -iNfft/2:iNfft/2-1;
            plot(kAxis,m(1)*kAxis+m(2))
            hold off
        end

        % Compensate Phase correction
        kAxis = -iNfft/2:iNfft/2-1;
        mPhaseEst = m(1)*kAxis+m(2);
        Rlk_out(l,:) = Rlk_in(l,:) .* exp(-1j*mPhaseEst);
    end
end
