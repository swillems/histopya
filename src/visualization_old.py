#!venv/bin/python


bottom_space = 0.3
button_count = 5
axcolor = 'lightgoldenrodyellow'


MIN_REP_COUNT_AX = plt.axes([0.65, bottom_space / button_count, 0.3, bottom_space / (button_count + 1)], facecolor=axcolor)
MAX_REP_COUNT_AX = plt.axes([0.65, 2 * bottom_space / button_count, 0.3, bottom_space / (button_count + 1)], facecolor=axcolor)
LOG_FDR_AX = plt.axes([0.65, 3 * bottom_space / button_count, 0.3, bottom_space / (button_count + 1)], facecolor=axcolor)

MIN_REP_COUNT_SLIDER = Slider(MIN_REP_COUNT_AX, "MIN_REP_COUNT", parameters["SIGNAL_COUNT_THRESHOLD"], parameters["SAMPLE_COUNT"], MIN_REP_COUNT, valstep=1)
MAX_REP_COUNT_SLIDER = Slider(MAX_REP_COUNT_AX, "MAX_REP_COUNT", parameters["SIGNAL_COUNT_THRESHOLD"], parameters["SAMPLE_COUNT"], MAX_REP_COUNT, valstep=1)
LOG_FDR_SLIDER = Slider(LOG_FDR_AX, "LOG_FDR", -5, 0, valinit=LOG_FDR, valstep=0.01, valfmt="%1.3f")


ax_network = plt.axes(
    [
        0.45,
        bottom_space / button_count,
        0.1,
        (2 * bottom_space / button_count + bottom_space / (button_count + 1)) / 2
    ],
    facecolor=axcolor
)
network_radio = RadioButtons(ax_network, ('None', 'Selected', 'All'), active=0)


def update_network_plot(label):
    global NETWORK_VISIBLE
    NETWORK_VISIBLE = label
    plt_bg_nx()


network_radio.on_clicked(update_network_plot)

ax_label_selection = plt.axes(
    [
        0.45,
        bottom_space / button_count + (2 * bottom_space / button_count + bottom_space / (button_count + 1)) / 2,
        0.1,
        (2 * bottom_space / button_count + bottom_space / (button_count + 1)) / 2
    ],
    facecolor=axcolor
)
label_selection_radio = RadioButtons(ax_label_selection, ('None', 'm/z', 'Peptide'), active=0)


def update_label_selection_plot(label):
    global LABEL_SELECTION
    LABEL_SELECTION = label
    plt_anchor_selection()
    plt_intensities()


label_selection_radio.on_clicked(update_label_selection_plot)

ax_expand_neighbors = plt.axes([0.3, 3 * bottom_space / button_count, 0.1, bottom_space / (button_count + 1)], facecolor=axcolor)
neighbor_button = Button(ax_expand_neighbors, 'Expand neighbors', color=axcolor, hovercolor='0.975')


def expand_neighbors(val):
    selection = browser.anchors["ION_COUNT"] >= MIN_REP_COUNT
    selection &= browser.anchors["ION_COUNT"] <= MAX_REP_COUNT
    dt_low, dt_high = plt.ylim()
    rt_low, rt_high = plt.xlim()
    selection &= browser.anchors["DT"] <= dt_high
    selection &= browser.anchors["DT"] >= dt_low
    selection &= browser.anchors["RT"] <= rt_high
    selection &= browser.anchors["RT"] >= rt_low
    selection = np.nonzero(selection)
    n = np.unique(browser.neighbors[self.SELECTED_ANCHORS].indices)
    n = n[np.isin(n, selection)]
    for anchor_ind in n:
        if anchor_ind not in self.SELECTED_ANCHORS:
            self.SELECTED_ANCHORS.append(anchor_ind)
    # print("Peptides: ", np.unique(anchor_peptide_match_counts[self.SELECTED_ANCHORS].indices, return_counts=True))
    plt_bg_nx()
    plt_anchor_selection()
    plt_intensities()


neighbor_button.on_clicked(expand_neighbors)

ax_refresh = plt.axes([0.3, bottom_space / button_count, 0.1, bottom_space / (button_count + 1)], facecolor=axcolor)
refresh_button = Button(ax_refresh, 'Refresh', color=axcolor, hovercolor='0.975')


def refresh(val):
    global self.SELECTED_ANCHORS
    selection = browser.anchors["ION_COUNT"] >= MIN_REP_COUNT
    selection &= browser.anchors["ION_COUNT"] <= MAX_REP_COUNT
    dt_low, dt_high = plt.ylim()
    rt_low, rt_high = plt.xlim()
    selection &= browser.anchors["DT"] <= dt_high
    selection &= browser.anchors["DT"] >= dt_low
    selection &= browser.anchors["RT"] <= rt_high
    selection &= browser.anchors["RT"] >= rt_low
    selection = np.nonzero(selection)
    self.SELECTED_ANCHORS = [self.SELECTED_ANCHORS[i] for i in np.flatnonzero(np.isin(self.SELECTED_ANCHORS, selection))]
    plt_bg_nx()
    plt_bg()
    plt_fg()
    plt_anchor_selection()
    plt_intensities()


refresh_button.on_clicked(refresh)

ax_clear = plt.axes([0.3, 2 * bottom_space / button_count, 0.1, bottom_space / (button_count + 1)], facecolor=axcolor)
clear_button = Button(ax_clear, 'Clear selection', color=axcolor, hovercolor='0.975')


def clear_selection(val):
    del self.SELECTED_ANCHORS[:]
    plt_bg_nx()
    plt_anchor_selection()ection)



def anchors_update(val):
    global MIN_REP_COUNT
    global MAX_REP_COUNT
    global LOG_FDR
    MIN_REP_COUNT = MIN_REP_COUNT_SLIDER.val
    MAX_REP_COUNT = MAX_REP_COUNT_SLIDER.val
    LOG_FDR = LOG_FDR_SLIDER.val
    refresh(None)


MIN_REP_COUNT_SLIDER.on_changed(anchors_update)
MAX_REP_COUNT_SLIDER.on_changed(anchors_update)
LOG_FDR_SLIDER.on_changed(anchors_update)


def onclick(event):
    if not event.dblclick:
        return
    rt = event.xdata
    dt = event.ydata
    selection = browser.anchors["ION_COUNT"] >= MIN_REP_COUNT
    selection &= browser.anchors["ION_COUNT"] <= MAX_REP_COUNT
    dt_low, dt_high = plt.ylim()
    rt_low, rt_high = plt.xlim()
    selection &= browser.anchors["DT"] <= dt_high
    selection &= browser.anchors["DT"] >= dt_low
    selection &= browser.anchors["RT"] <= rt_high
    selection &= browser.anchors["RT"] >= rt_low
    rt_distance = (browser.anchors["RT"][selection] - rt) / (rt_high - rt_low)
    dt_distance = (browser.anchors["DT"][selection] - dt) / (dt_high - dt_low)
    distance = np.sqrt(rt_distance**2 + dt_distance**2)
    min_distance_ind = np.argmin(distance)
    anchor_ind = np.flatnonzero(selection)[min_distance_ind]
    if event.key != "control":
        del self.SELECTED_ANCHORS[:]
    if anchor_ind in self.SELECTED_ANCHORS:
        sel_i = self.SELECTED_ANCHORS.index(anchor_ind)
        del self.SELECTED_ANCHORS[sel_i]
    else:
        self.SELECTED_ANCHORS.append(anchor_ind)
    # print("Peptides: ", np.unique(anchor_peptide_match_counts[self.SELECTED_ANCHORS].indices, return_counts=True))
    plt_bg_nx()
    plt_anchor_selection()
    plt_intensities()


plt_intensities = browser.plot_intensities()

ind = fig.canvas.mpl_connect('button_press_event', onclick)

    plt_intensities()


clear_button.on_clicked(clear_selection)



def anchors_update(val):
    global MIN_REP_COUNT
    global MAX_REP_COUNT
    global LOG_FDR
    MIN_REP_COUNT = MIN_REP_COUNT_SLIDER.val
    MAX_REP_COUNT = MAX_REP_COUNT_SLIDER.val
    LOG_FDR = LOG_FDR_SLIDER.val
    refresh(None)


MIN_REP_COUNT_SLIDER.on_changed(anchors_update)
MAX_REP_COUNT_SLIDER.on_changed(anchors_update)
LOG_FDR_SLIDER.on_changed(anchors_update)


def onclick(event):
    if not event.dblclick:
        return
    rt = event.xdata
    dt = event.ydata
    selection = browser.anchors["ION_COUNT"] >= MIN_REP_COUNT
    selection &= browser.anchors["ION_COUNT"] <= MAX_REP_COUNT
    dt_low, dt_high = plt.ylim()
    rt_low, rt_high = plt.xlim()
    selection &= browser.anchors["DT"] <= dt_high
    selection &= browser.anchors["DT"] >= dt_low
    selection &= browser.anchors["RT"] <= rt_high
    selection &= browser.anchors["RT"] >= rt_low
    rt_distance = (browser.anchors["RT"][selection] - rt) / (rt_high - rt_low)
    dt_distance = (browser.anchors["DT"][selection] - dt) / (dt_high - dt_low)
    distance = np.sqrt(rt_distance**2 + dt_distance**2)
    min_distance_ind = np.argmin(distance)
    anchor_ind = np.flatnonzero(selection)[min_distance_ind]
    if event.key != "control":
        del self.SELECTED_ANCHORS[:]
    if anchor_ind in self.SELECTED_ANCHORS:
        sel_i = self.SELECTED_ANCHORS.index(anchor_ind)
        del self.SELECTED_ANCHORS[sel_i]
    else:
        self.SELECTED_ANCHORS.append(anchor_ind)
    # print("Peptides: ", np.unique(anchor_peptide_match_counts[self.SELECTED_ANCHORS].indices, return_counts=True))
    plt_bg_nx()
    plt_anchor_selection()
    plt_intensities()


plt_intensities = browser.plot_intensities()

ind = fig.canvas.mpl_connect('button_press_event', onclick)
