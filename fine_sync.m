function [Rlk_out] = fine_sync(Rlk_in, Plk, iNfft, iNOfSymbolsPerFrame, SwitchDemoSync)
    % FINE_SYNC Linear phase correction using pilot tones
    % Performs fine synchronization by estimating and correcting linear phase errors
    %
    % Inputs:
    %   Rlk_in - Received OFDM frame in frequency domain
    %   Plk - Pilot pattern matrix
    %   iNfft - FFT length
    %   iNOfSymbolsPerFrame - Number of symbols per frame
    %   SwitchDemoSync - Enable debug plots
    %
    % Outputs:
    %   Rlk_out - Phase-corrected OFDM frame

    % Copy input to output
    Rlk_out = Rlk_in;

    % Calculate number of frames
    iNOfFrames = floor(size(Rlk_in, 1) / iNOfSymbolsPerFrame);

    % For each frame
    for i = 1:iNOfFrames
        % Get current frame's symbols
        iFrameStart = (i-1) * iNOfSymbolsPerFrame + 1;
        RlkTemp = Rlk_in(iFrameStart:iFrameStart+iNOfSymbolsPerFrame-1, :);

        % Calculate phase difference between received signal and known pilots
        mPhase = angle(conj(RlkTemp) .* Plk);

        % Plot phase information if debug enabled
        if SwitchDemoSync
            figure(301)
            mesh([-iNfft/2:iNfft/2-1], [1:size(Plk,1)], unwrap(mPhase))
            ylabel('Symbol l')
            xlabel('Subchannel k')
            zlabel('Phase \Phi(l,k)')
        end

        % For each symbol in the frame
        for l = 1:size(Plk, 1)
            % Plot per-symbol phase if debug enabled
            if SwitchDemoSync
                figure(302)
                stem([-iNfft/2:iNfft/2-1], unwrap(mPhase(l,:)))
                xlabel('Subchannel k')
                ylabel(['Phase \Phi(' num2str(l) ',k)'])
            end

            % Linear regression of phases
            kPilots = [find(Plk(l,:) ~= 0) - iNfft/2 - 1].';
            V = [kPilots ones(size(kPilots))];
            vPhi = unwrap(mPhase(l, Plk(l,:) ~= 0)).';

            % Solve linear system for phase offset parameters
            m = V\vPhi;

            % Create subchannel axis
            kAxis = -iNfft/2:iNfft/2-1;

            % Plot regression line if debug enabled
            if SwitchDemoSync
                hold on
                plot(kAxis, m(1)*kAxis+m(2))
                hold off
            end

            % Calculate phase correction for each subchannel
            mPhaseEst = m(1)*kAxis + m(2);

            % Apply phase correction to current symbol
            Rlk_out(iFrameStart+l-1, :) = Rlk_out(iFrameStart+l-1, :) .* exp(-1*j*mPhaseEst);
        end
    end
end
