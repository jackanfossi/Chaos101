%  TRNG ANALYSIS

fname = input('Name of the file .bin (ex: chua_data.bin): ', 's');

% Same folder
fpath = fileparts(mfilename('fullpath'));

fid = fopen(fullfile(fpath, fname), 'rb');
if fid == -1
    fprintf('ERRORE: file "%s" non trovato nella cartella dello script.\n', fname);
    return;
end
raw_bytes = fread(fid, Inf, 'uint8');
fclose(fid);
bits = de2bi(raw_bytes, 8, 'left-msb');
bits = bits(:)';

N_bytes = length(raw_bytes);
N_bits  = length(bits);

%  1. ENTROPY
p1    = mean(bits);
p0    = 1 - p1;
H_bit = -p1*log2(p1+eps) - p0*log2(p0+eps);

counts_byte = histcounts(raw_bytes, 0:256);
p_byte      = counts_byte / N_bytes;
H_byte      = -sum(p_byte .* log2(p_byte + eps));

fprintf('=== ENTROPIA ===\n');
fprintf('Bit  per bit : %.6f / 1.000000\n', H_bit);
fprintf('Byte per byte: %.6f / 8.000000\n\n', H_byte);

%  2. CHI-QUADRO
expected_bit  = N_bits / 2;
chi2_bit      = ((sum(bits==0) - expected_bit)^2 + ...
                 (sum(bits==1) - expected_bit)^2) / expected_bit;
p_chi2_bit    = 1 - chi2cdf(chi2_bit, 1);

expected_byte = N_bytes / 256;
chi2_byte     = sum((counts_byte - expected_byte).^2 / expected_byte);
p_chi2_byte   = 1 - chi2cdf(chi2_byte, 255);

fprintf('=== CHI-QUADRO ===\n');
fprintf('Bit  : chi2=%.4f  p=%.2f%%\n', chi2_bit,  p_chi2_bit*100);
fprintf('Byte : chi2=%.4f  p=%.2f%%\n\n', chi2_byte, p_chi2_byte*100);

%  3. MEDIA ARITMETICA
mean_bit  = mean(bits);
mean_byte = mean(raw_bytes);

fprintf('=== MEDIA ARITMETICA ===\n');
fprintf('Bit  : %.4f  (ideale = 0.5)\n',    mean_bit);
fprintf('Byte : %.4f  (ideale = 127.5)\n\n', mean_byte);

%  4. MONTE CARLO PI
n_groups = floor(N_bytes / 6);
idx      = reshape(raw_bytes(1:n_groups*6), 6, n_groups)';
X        = idx(:,1).*65536 + idx(:,2).*256 + idx(:,3);
Y        = idx(:,4).*65536 + idx(:,5).*256 + idx(:,6);
max_val  = 2^24 - 1;
dist     = sqrt((X/max_val).^2 + (Y/max_val).^2);
hits     = sum(dist <= 1);
pi_mc    = 4 * hits / n_groups;
pi_err   = abs(pi_mc - pi) / pi * 100;

fprintf('=== MONTE CARLO PI ===\n');
fprintf('Pi stimato : %.9f\n', pi_mc);
fprintf('Pi reale   : %.9f\n', pi);
fprintf('Errore     : %.2f%%\n\n', pi_err);

%  5. CORRELAZIONE SERIALE
scf_bit  = corr(bits(1:end-1)', bits(2:end)');
scf_byte = corr(double(raw_bytes(1:end-1)), double(raw_bytes(2:end)));

fprintf('=== CORRELAZIONE SERIALE ===\n');
fprintf('Bit  : %.6f  (ideale = 0.0)\n', scf_bit);
fprintf('Byte : %.6f  (ideale = 0.0)\n\n', scf_byte);


dark  = [1.00 1.00 1.00];
panel = [0.97 0.97 0.97];
cyan  = [0.00 0.45 0.74];
purp  = [0.49 0.18 0.56];
yell  = [0.85 0.55 0.00];
red_c = [0.85 0.15 0.15];
grn   = [0.10 0.65 0.20];
white = [0.10 0.10 0.10];
muted = [0.30 0.30 0.30];


%  Distribution Bit plot
figure('Name','1 - Distribuzione Bit','Color',dark,'Position',[100 500 500 380]);
ax = gca; style_ax(ax, panel, white, muted);
bar([0 1],[sum(bits==0) sum(bits==1)],0.5,'FaceColor',purp,'EdgeColor','none');
hold on;
yline(N_bits/2,'--','Color',cyan,'LineWidth',1.5);
title('Distribuzione Bit'); xlabel('Valore bit'); ylabel('Conteggio');
xticklabels({'0','1'});
legend({'Conteggi','Ideale'},'TextColor',white,'Color',panel,'EdgeColor',muted);
grid on;

% Distribuzione Byte plot
figure('Name','2 - Distribuzione Byte','Color',dark,'Position',[620 500 500 380]);
ax = gca; style_ax(ax, panel, white, muted);
bar(0:255,counts_byte,1,'FaceColor',cyan,'EdgeColor','none','FaceAlpha',0.85);
hold on;
yline(expected_byte,'--','Color',yell,'LineWidth',1.5);
title('Distribuzione Byte (0-255)'); xlabel('Valore byte'); ylabel('Conteggio');
legend({'Conteggi','Ideale'},'TextColor',white,'Color',panel,'EdgeColor',muted);
grid on;


%  Entropy plot
figure('Name','3 - Entropia','Color',dark,'Position',[100 60 500 380]);
ax = gca; style_ax(ax, panel, white, muted);
vals = [H_bit, H_byte/8];
b = bar(vals,0.45,'FaceColor','flat','EdgeColor','none');
b.CData(1,:) = purp; b.CData(2,:) = cyan;
hold on;
yline(1.0,'--','Color',yell,'LineWidth',1.5);
set(gca,'XTickLabel',{'Bit/bit','Byte/byte (÷8)'},'YLim',[0.98 1.002]);
title('Entropia normalizzata'); ylabel('Entropia / max');
text(1, H_bit+0.0003,    sprintf('%.6f',H_bit),   'Color',white,'HorizontalAlignment','center','FontSize',9);
text(2, H_byte/8+0.0003, sprintf('%.6f',H_byte/8),'Color',white,'HorizontalAlignment','center','FontSize',9);
grid on;

%  Entropy
%figure('Name','4 - Entropia nel tempo','Color',dark,'Position',[620 60 500 380]);
%ax = gca; style_ax(ax, panel, white, muted);
%window = 1000; step = 200;
%%idxs  = window:step:N_bits;
% = arrayfun(@(i) ...
 %   -mean(bits(i-window+1:i))*log2(mean(bits(i-window+1:i))+eps) ...
 %   -(1-mean(bits(i-window+1:i)))*log2(1-mean(bits(i-window+1:i))+eps), idxs);
%plot(idxs, H_run,'Color',cyan,'LineWidth',1.2);
%hold on;
%yline(1.0,'--','Color',yell,'LineWidth',1.2);
%title('Entropia locale nel tempo'); xlabel('Posizione bit'); ylabel('H [bit/bit]');
%ylim([0.95 1.005]); grid on;

%  Monte Carlo Pi scatter
figure('Name','5 - Monte Carlo Pi','Color',dark,'Position',[1140 500 500 500]);
ax = gca; style_ax(ax, panel, white, muted);

% Coordinates [-1, +1] for the entire circle
n_show = min(5000, n_groups);
Xs = (X(1:n_show)/max_val)*2 - 1;
Ys = (Y(1:n_show)/max_val)*2 - 1;
d_show = sqrt(Xs.^2 + Ys.^2);
in_c = d_show<=1; out_c = ~in_c;

scatter(Xs(in_c),  Ys(in_c),  1, repmat(cyan, sum(in_c), 1),  'filled'); hold on;
scatter(Xs(out_c), Ys(out_c), 1, repmat(red_c,sum(out_c),1), 'filled');

% Entire circle
th = linspace(0, 2*pi, 500);
plot(cos(th), sin(th), '--', 'Color', yell, 'LineWidth', 2);

% Square
plot([-1 1 1 -1 -1], [-1 -1 1 1 -1], '-', 'Color', muted, 'LineWidth', 1.5);

title(sprintf('Monte Carlo  pi=%.6f  err=%.2f%%', pi_mc, pi_err));
xlabel('X'); ylabel('Y');
axis equal; axis([-1.05 1.05 -1.05 1.05]); grid on;
legend({'Dentro','Fuori','Cerchio'}, 'TextColor',white,'Color',panel,'EdgeColor',muted);

%  Autocorrelazione Bit
figure('Name','6 - Autocorrelazione','Color',dark,'Position',[1140 60 500 380]);
ax = gca; style_ax(ax, panel, white, muted);
max_lag = 100;
n_ac    = min(50000, N_bits);
[ac, lags] = xcorr(double(bits(1:n_ac))-mean_bit, max_lag, 'normalized');
plot(lags, ac,'Color',purp,'LineWidth',1); hold on;
conf = 1.96/sqrt(n_ac);
yline( conf,'--','Color',yell,'LineWidth',1);
yline(-conf,'--','Color',yell,'LineWidth',1);
yline(0,'Color',white,'LineWidth',0.5);
title('Autocorrelazione Bit'); xlabel('Lag'); ylabel('Correlazione');
legend({'ACF','Intervallo 95%'},'TextColor',white,'Color',panel,'EdgeColor',muted);
grid on;

% =========================================================================
%  FIGURA 7 — Scatter byte consecutivi
% =========================================================================
figure('Name','7 - Scatter Byte','Color',dark,'Position',[100 500 500 380]);
ax = gca; style_ax(ax, panel, white, muted);
n_sc = min(5000, N_bytes-1);
scatter(double(raw_bytes(1:n_sc)), double(raw_bytes(2:n_sc+1)), ...
    1, repmat(cyan,n_sc,1),'filled','MarkerFaceAlpha',0.3);
title(sprintf('Byte[n] vs Byte[n+1]  SCF=%.4f',scf_byte));
xlabel('Byte[n]'); ylabel('Byte[n+1]');
axis([0 255 0 255]); grid on;

% =========================================================================
%  FIGURA 8 — Convergenza Monte Carlo Pi
% =========================================================================
figure('Name','5 - Monte Carlo Pi','Color',[1 1 1],'Position',[1140 500 560 520]);
ax = gca;
set(ax, 'Color', [1 1 1], 'XColor', [0.2 0.2 0.2], 'YColor', [0.2 0.2 0.2], ...
    'FontName','Arial', 'FontSize', 10, 'GridColor',[0.8 0.8 0.8], 'GridAlpha',0.5);
ax.Title.Color = [0.1 0.1 0.1];
ax.XLabel.Color = [0.2 0.2 0.2];
ax.YLabel.Color = [0.2 0.2 0.2];

% Coordinate in [-1, +1]
n_show = min(10000, n_groups);
Xs = (X(1:n_show)/max_val)*2 - 1;
Ys = (Y(1:n_show)/max_val)*2 - 1;
d_show = sqrt(Xs.^2 + Ys.^2);
in_c  = d_show <= 1;
out_c = ~in_c;

% Punti fuori (rossi) prima, dentro (blu) sopra
scatter(Xs(out_c), Ys(out_c), 4, [0.85 0.15 0.15], 'filled', ...
    'MarkerFaceAlpha', 0.6); hold on;
scatter(Xs(in_c),  Ys(in_c),  4, [0.15 0.35 0.85], 'filled', ...
    'MarkerFaceAlpha', 0.6);

% Cerchio intero
th = linspace(0, 2*pi, 600);
plot(cos(th), sin(th), 'k-', 'LineWidth', 1.8);

% Quadrato
plot([-1 1 1 -1 -1], [-1 -1 1 1 -1], '-', ...
    'Color',[0.4 0.4 0.4], 'LineWidth', 1.2);

title(sprintf('Monte Carlo Simulation for \\pi'), ...
    'FontSize', 13, 'FontWeight','bold', 'FontName','Arial');
text(-0.95, 1.08, sprintf('Estimated value of \\pi with %d samples: %.4f', ...
    n_show, pi_mc), 'FontSize', 9, 'FontName','Arial', 'Color',[0.2 0.2 0.2]);
xlabel('X'); ylabel('Y');
axis equal; axis([-1.05 1.05 -1.05 1.05]); grid on;
legend({'Fuori','Dentro','Cerchio'}, 'Location','southeast', ...
    'FontSize',9, 'FontName','Arial');


% =========================================================================
%  Bitmap visuale della sequenza
figure('Name','Bitmap TRNG','Color',dark,'Position',[100 100 500 500]);
ax = gca; style_ax(ax, panel, white, muted);

% Calcola dimensioni quadrate
lato = floor(sqrt(N_bits));
bits_quadro = bits(1:lato*lato);          % prendi solo i bit necessari
bitmap = reshape(bits_quadro, lato, lato); % matrice quadrata

imagesc(bitmap);
colormap(gray);
axis equal; axis tight; axis off;
title(sprintf('Bitmap TRNG  %d x %d pixel', lato, lato), ...
      'Color', white, 'FontName', 'Courier New', 'FontSize', 10);


%  FIPS 140-2 TEST
%  At least 20000 bit
fprintf('=== FIPS 140-2 ===\n');

if N_bits < 20000
    fprintf('ERRORE: servono almeno 20000 bit (hai %d)\n', N_bits);
else
    fips_bits = bits(1:20000);

    % ── Monobit Test ──────────────────────────────────────────────────
    % Conta gli 1 nella sequenza — deve essere tra 9725 e 10275
    monobit = sum(fips_bits);
    monobit_pass = monobit >= 9725 && monobit <= 10275;
    fprintf('Monobit    : %d uni  →  %s\n', monobit, ...
        string_pass(monobit_pass));

    % ── Poker Test ────────────────────────────────────────────────────
    % Divide i 20000 bit in gruppi da 4, conta le 16 combinazioni possibili
    % La statistica deve essere tra 2.16 e 46.17
    poker_groups = reshape(fips_bits(1:20000), 4, 5000)';
    poker_vals   = poker_groups * [8;4;2;1];   % converti ogni gruppo in 0-15
    poker_counts = histcounts(poker_vals, -0.5:1:15.5);
    poker_stat   = (16/5000) * sum(poker_counts.^2) - 5000;
    poker_pass   = poker_stat >= 2.16 && poker_stat <= 46.17;
    fprintf('Poker      : %.4f  →  %s\n', poker_stat, ...
        string_pass(poker_pass));

    % ── Runs Test ─────────────────────────────────────────────────────
    % Conta le sequenze consecutive di 0 o 1 (run)
    % Ogni lunghezza deve rientrare in un range specifico
    runs_limits = [2315 2685; 1114 1386; 527 723; 240 384; 103 209; 103 209];
    runs_counts = zeros(1,6);
    current_len = 1;
    for i = 2:20000
        if fips_bits(i) == fips_bits(i-1)
            current_len = current_len + 1;
        else
            idx = min(current_len, 6);
            runs_counts(idx) = runs_counts(idx) + 1;
            current_len = 1;
        end
    end
    % ultimo run
    idx = min(current_len, 6);
    runs_counts(idx) = runs_counts(idx) + 1;

    runs_pass = true;
    for r = 1:6
        ok = runs_counts(r) >= runs_limits(r,1) && ...
             runs_counts(r) <= runs_limits(r,2);
        if ~ok, runs_pass = false; end
    end
    fprintf('Runs       : %s\n', string_pass(runs_pass));
    for r = 1:6
        ok = runs_counts(r) >= runs_limits(r,1) && ...
             runs_counts(r) <= runs_limits(r,2);
        fprintf('  Lung.%d : %4d  (range %d-%d)  %s\n', r, ...
            runs_counts(r), runs_limits(r,1), runs_limits(r,2), ...
            string_pass(ok));
    end

    % ── Long Run Test ─────────────────────────────────────────────────
    % Non deve esserci nessuna sequenza consecutiva >= 26 dello stesso bit
    max_run = 1; cur = 1;
    for i = 2:20000
        if fips_bits(i) == fips_bits(i-1)
            cur = cur + 1;
            if cur > max_run, max_run = cur; end
        else
            cur = 1;
        end
    end
    longrun_pass = max_run < 26;
    fprintf('Long Run   : max=%d  →  %s\n\n', max_run, ...
        string_pass(longrun_pass));

    % ── Risultato finale ──────────────────────────────────────────────
    fips_pass = monobit_pass && poker_pass && runs_pass && longrun_pass;
    if fips_pass
        fprintf('FIPS 140-2 : PASSED ✅\n\n');
    else
        fprintf('FIPS 140-2 : FAILED ❌\n\n');
    end

    % ── Grafico FIPS ──────────────────────────────────────────────────
    figure('Name','FIPS 140-2','Color',dark,'Position',[200 200 500 420]);
    ax = gca; style_ax(ax, panel, white, muted);
    axis off;

    tests  = {'Monobit','Poker','Runs','Long Run','FIPS 140-2 totale'};
    passed = [monobit_pass, poker_pass, runs_pass, longrun_pass, fips_pass];
    values = {sprintf('%d / [9725-10275]', monobit), ...
              sprintf('%.2f / [2.16-46.17]', poker_stat), ...
              sprintf('%s', string_pass(runs_pass)), ...
              sprintf('max run = %d / <26', max_run), ...
              ''};

    title(ax, 'FIPS 140-2  (20.000 bit)', 'Color',white, ...
          'FontName','Courier New','FontSize',11,'FontWeight','bold');

    for i = 1:5
        col = grn*double(passed(i)) + red_c*double(~passed(i));
        text(0.02, 1-i*0.17, sprintf('●  %s', tests{i}), ...
             'Color',col,'FontName','Courier New','FontSize',10,...
             'Units','normalized','Parent',ax,'FontWeight','bold');
        text(0.98, 1-i*0.17, values{i}, ...
             'Color',white,'FontName','Courier New','FontSize',9,...
             'HorizontalAlignment','right','Units','normalized','Parent',ax);
    end
end


function style_ax(ax, panel, white, muted)
    set(ax, 'Color', panel, 'XColor', muted, 'YColor', muted, ...
    'GridColor', [0.8 0.8 0.8], 'GridAlpha', 0.6, ... 
        'TitleHorizontalAlignment','left');
    ax.Title.Color  = white;
    ax.XLabel.Color = muted;
    ax.YLabel.Color = muted;
end

function s = string_pass(ok)
    if ok
        s = 'PASS';
    else
        s = 'FAIL';
    end
end