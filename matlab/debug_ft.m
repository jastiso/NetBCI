 cfg = [];
               if strcmp('mag',sens)
                cfg.layout = [top_dir, '/layouts/neuromag306', sens, '.lay'];
               else
                   cfg.layout = [top_dir, '/layouts/neuromag306planar.lay'];
               end
                plot_data.powsctrm = sg_deg;
               plot_data.labels = labels;
               plot_data.dimord = 'chan_freq';
               plot_data.freq = 6;
               plot_data.cfg = [];
               figure; 
               ft_topoplotER(cfg,plot_data); colorbar