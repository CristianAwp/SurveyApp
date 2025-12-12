// --- MATH UTILITIES ---
const MathUtils = {
    // Convert Grads to Radians (1 circle = 400g = 2PI rad)
    gradsToRad: (g) => g * (Math.PI / 200),
    
    // Convert Radians to Grads
    radToGrads: (r) => r * (200 / Math.PI),

    // Normalize angle to 0-400 grads
    normalizeGrads: (g) => {
        g = g % 400;
        if (g < 0) g += 400;
        return g;
    },

    // Calculate Horizontal Distance from Stadia
    // D = K * (Top - Bottom) * sin^2(Zenith)
    calcDist: (top, bot, vInput, k) => {
        const intercept = top - bot;
        
        let z = vInput;
        if(State.settings.vMode === 'alpha') {
            z = 100 - vInput;
        }
        
        const zRad = MathUtils.gradsToRad(z);
        return k * intercept * Math.pow(Math.sin(zRad), 2);
    },

    // Calculate Vertical Difference (dH) from Station Axis to Target Axis (at mid wire)
    // V = 0.5 * K * (top - bot) * sin(2*Zenith)
    // Note: This dH is from Instrument Axis to Target axis (where mid wire is).
    // Total Height Diff = dH + InstHeight - MidWireReading
    calcVertDiffComponent: (top, bot, vInput, k) => {
        const intercept = top - bot;
        
        let z = vInput;
        if(State.settings.vMode === 'alpha') {
            z = 100 - vInput;
        }

        const zRad = MathUtils.gradsToRad(z);
        return 0.5 * k * intercept * Math.sin(2 * zRad);
    },

    // Coordinate Calculation
    polar: (e, n, azGrads, dist) => {
        const azRad = MathUtils.gradsToRad(azGrads);
        return {
            e: e + dist * Math.sin(azRad),
            n: n + dist * Math.cos(azRad)
        };
    },

    // Azimuth between two coordinates
    azimuth: (e1, n1, e2, n2) => {
        const de = e2 - e1;
        const dn = n2 - n1;
        let azRad = Math.atan2(de, dn);
        if (azRad < 0) azRad += 2 * Math.PI;
        return MathUtils.radToGrads(azRad);
    },

    dist2D: (e1, n1, e2, n2) => {
        return Math.sqrt(Math.pow(e2-e1, 2) + Math.pow(n2-n1, 2));
    },

    // Shoelace formula for polygon area
    calcPolygonArea: (points) => {
        let area = 0;
        const j = points.length - 1;
        for (let i = 0; i < points.length; i++) {
            const next = (i + 1) % points.length;
            area += points[i].e * points[next].n;
            area -= points[i].n * points[next].e;
        }
        return Math.abs(area) / 2;
    },

    calcPerimeter: (points) => {
        let p = 0;
        for (let i = 0; i < points.length; i++) {
            const next = (i + 1) % points.length;
            p += MathUtils.dist2D(points[i].e, points[i].n, points[next].e, points[next].n);
        }
        return p;
    },

    // --- RESECTION (LEAST SQUARES) ---
    // Solver for 3 unknowns (dE, dN, dOri) using observations to multiple points
    solveResection: (obs) => {
        // obs: Array of { e, n, dist, hz } (Target Coords, Measured Dist, Measured Hz in Grads)
        // Initial approximation: Use the first point. 
        // Assume Stn = Pt1 coords (bad guess but iterative logic handles it if robust, 
        // better: use first point and arbitrary orientation)
        
        // Initial Guess:
        // Stn E, N = Pt1 E, N (offset by distance). 
        // Let's assume Stn is at (Pt1.E - dist, Pt1.N) for start.
        // Or simpler: Just start at Average of target coords? No, that might be inside the polygon.
        // Best simple start: "Polar" from Pt 1 reversed.
        
        let stnE = obs[0].e;
        let stnN = obs[0].n; // Very rough start (station at target 1)
        // Better start: If we have dist to pt1, assume we are 'dist' away to South.
        stnN -= obs[0].dist; 
        
        let ori = 0; // Initial orientation guess
        
        // Iterations
        for(let iter=0; iter<10; iter++) {
            let A = []; // Jacobian
            let L = []; // Residuals (Observed - Computed)
            
            // For each observation
            for(let i=0; i<obs.length; i++) {
                const pt = obs[i];
                const dx = pt.e - stnE;
                const dy = pt.n - stnN;
                const distCalc = Math.sqrt(dx*dx + dy*dy);
                const azCalcRad = Math.atan2(dx, dy); // North Azimuth in Radians
                let azCalc = azCalcRad * (200/Math.PI);
                if(azCalc < 0) azCalc += 400;

                // 1. Distance Equation
                // D_calc + (dD/dE)dE + (dD/dN)dN = D_obs
                // v = D_obs - D_calc
                // Coeffs: -sin(Az), -cos(Az)  (Derivatives w.r.t Station pos)
                const sinAz = Math.sin(azCalcRad);
                const cosAz = Math.cos(azCalcRad);
                
                A.push([ -sinAz, -cosAz, 0 ]); // 0 for dOri impact on distance
                L.push(pt.dist - distCalc);

                // 2. Angle Equation
                // Az_calc + dAz - Ori - dOri = Hz_obs
                // Hz_obs = Az_calc - Ori
                // Residual = (Hz_obs + Ori) - Az_calc
                // Actually: (Az_calc - Ori) - Hz_obs = v
                // Linearization:
                // dAz/dE * dE + dAz/dN * dN - 1 * dOri = Hz_obs - (Az_calc - Ori)
                // dAz/dE = cos(Az)/D, dAz/dN = -sin(Az)/D
                
                // Let's stick to standard: Observed = F(Params)
                // Hz_obs = Azimuth(Stn, Pt) - Orientation
                // Linearized:
                // Hz_obs = Az0 - Ori0 + (dAz/dE)dE + (dAz/dN)dN + (-1)dOri
                // L = Hz_obs - (Az0 - Ori0)
                // A row: [ cosAz/distCalc, -sinAz/distCalc, -1 ]
                
                // Handle 400g wrapping for residuals
                let angRes = pt.hz - (azCalc - ori);
                while(angRes > 200) angRes -= 400;
                while(angRes < -200) angRes += 400;
                
                // Convert angular residual to linear approximation (radians * dist)? 
                // No, standard LSA keeps units consistent if we weight them. 
                // For simplicity here, treat angle residual in Radians to match distance scale?
                // Or just solve in mixed units. 
                // Let's use Dist approx for weights or just simplistic unweighted.
                // Angle equation in Grads:
                
                // Note: To mix D (meters) and Angle (Grads), we need weights.
                // Simple approach: Multiply Angle row by Dist to make it "Transverse deviation in meters".
                // Eq * Dist:
                // Dist * (dAz/dE) = -cos(Az)  <-- Corrected Sign
                // Dist * (dAz/dN) = +sin(Az)  <-- Corrected Sign
                // Dist * (-1) * dOri = -Dist
                // RHS = Dist * AngleResidual(in Radians!)
                
                const angResRad = angRes * (Math.PI/200);
                
                // Weighted by Distance (Quasi-Metric):
                // dAz/dE = -cos(Az)/Dist, dAz/dN = sin(Az)/Dist
                // Multiplied by Dist:
                // [-cosAz, sinAz, -distCalc]
                
                A.push([ -cosAz, sinAz, -distCalc * (Math.PI/200) ]); 
                L.push(distCalc * angResRad); 
            }

            // Solve Normal Equations: (At A) x = (At L)
            const At = Matrix.transpose(A);
            const N = Matrix.multiply(At, A);
            const U = Matrix.multiplyVector(At, L);
            
            // Invert N (3x3)
            const Ninv = Matrix.invert3x3(N);
            if(!Ninv) return null; // Singular
            
            const x = Matrix.multiplyVector3x3(Ninv, U); // [dE, dN, dOri]
            
            stnE += x[0];
            stnN += x[1];
            ori += x[2];
            
            if(Math.abs(x[0]) < 0.001 && Math.abs(x[1]) < 0.001) break; // Converged
        }
        
        return { e: stnE, n: stnN, ori: MathUtils.normalizeGrads(ori) };
    }
};

const Matrix = {
    transpose: (A) => {
        return A[0].map((_, c) => A.map(r => r[c]));
    },
    multiply: (A, B) => {
        const result = new Array(A.length).fill(0).map(() => new Array(B[0].length).fill(0));
        return result.map((row, i) => {
            return row.map((val, j) => {
                return A[i].reduce((sum, elm, k) => sum + (elm * B[k][j]), 0)
            })
        })
    },
    multiplyVector: (A, v) => {
        return A.map(row => row.reduce((sum, val, i) => sum + val * v[i], 0));
    },
    multiplyVector3x3: (M, v) => {
        return [
            M[0][0]*v[0] + M[0][1]*v[1] + M[0][2]*v[2],
            M[1][0]*v[0] + M[1][1]*v[1] + M[1][2]*v[2],
            M[2][0]*v[0] + M[2][1]*v[1] + M[2][2]*v[2]
        ];
    },
    invert3x3: (m) => {
        const det = m[0][0] * (m[1][1] * m[2][2] - m[2][1] * m[1][2]) -
                    m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                    m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
        
        if (Math.abs(det) < 1e-8) return null;
        const invDet = 1 / det;
        
        return [
            [
                (m[1][1] * m[2][2] - m[2][1] * m[1][2]) * invDet,
                (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * invDet,
                (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * invDet
            ],
            [
                (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * invDet,
                (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * invDet,
                (m[1][0] * m[0][2] - m[0][0] * m[1][2]) * invDet
            ],
            [
                (m[1][0] * m[2][1] - m[2][0] * m[1][1]) * invDet,
                (m[2][0] * m[0][1] - m[0][0] * m[2][1]) * invDet,
                (m[0][0] * m[1][1] - m[1][0] * m[0][1]) * invDet
            ]
        ];
    }
};

// --- APP STATE ---
const State = {
    station: {
        e: 0, n: 0, z: 0, hi: 1.5,
        set: false
    },
    orientation: {
        azimuthCorrection: 0, // value to ADD to measured Hz to get Azimuth
        set: false
    },
    points: [], // Array of {id, e, n, z, code, ts}
    settings: {
        k: 100
    }
};

// --- DOM ELEMENTS ---
const UI = {
    tabs: document.querySelectorAll('.nav-btn'),
    contents: document.querySelectorAll('.tab-content'),
    
    // Status
    statusStn: document.getElementById('station-status'),
    statusOri: document.getElementById('orientation-status'),

    // Inputs - Setup
    kFactor: document.getElementById('k-factor'),
    distModeInputs: document.querySelectorAll('input[name="dist-mode"]'),
    vModeInputs: document.querySelectorAll('input[name="v-mode"]'),
    // Mode Switch
    modeKnown: document.getElementById('mode-known'),
    modeResection: document.getElementById('mode-resection'),
    setupKnown: document.getElementById('setup-known-point'),
    setupResection: document.getElementById('setup-resection'),

    // Known Point
    stnE: document.getElementById('stn-e'),
    stnN: document.getElementById('stn-n'),
    stnZ: document.getElementById('stn-z'),
    stnHi: document.getElementById('stn-hi'),
    bsE: document.getElementById('bs-e'),
    bsN: document.getElementById('bs-n'),
    bsHz: document.getElementById('bs-hz'),
    btnSetOri: document.getElementById('btn-set-orientation'),

    // Resection
    resectionList: document.getElementById('resection-list'),
    btnAddResTarget: document.getElementById('btn-add-res-target'),
    btnCalcRes: document.getElementById('btn-calc-resection'),
    resHi: document.getElementById('res-hi'),
    resResult: document.getElementById('resection-result'),
    resStnE: document.getElementById('res-stn-e'),
    resStnN: document.getElementById('res-stn-n'),
    resOri: document.getElementById('res-ori'),
    btnApplyRes: document.getElementById('btn-apply-resection'),

    // Inputs - Survey
    measPid: document.getElementById('meas-pid'),
    surveyStadiaDiv: document.getElementById('survey-stadia-inputs'),
    surveyManualDiv: document.getElementById('survey-manual-inputs'),
    measTop: document.getElementById('meas-top'),
    measMid: document.getElementById('meas-mid'),
    measBot: document.getElementById('meas-bot'),
    measSDist: document.getElementById('meas-sdist'),
    measHr: document.getElementById('meas-hr'),
    measHz: document.getElementById('meas-hz'),
    measV: document.getElementById('meas-v'),
    btnMeasure: document.getElementById('btn-measure'),
    resCard: document.getElementById('survey-result'),
    resDist: document.getElementById('res-dist'),
    resE: document.getElementById('res-e'),
    resN: document.getElementById('res-n'),
    resH: document.getElementById('res-h'),

    // Stakeout
    soTargetSelect: document.getElementById('so-target-select'),
    soTargetInfo: document.getElementById('so-target-info'),
    soStadiaDiv: document.getElementById('so-stadia-inputs'),
    soManualDiv: document.getElementById('so-manual-inputs'),
    soTgtE: document.getElementById('so-tgt-e'),
    soTgtN: document.getElementById('so-tgt-n'),
    soTop: document.getElementById('so-top'),
    soMid: document.getElementById('so-mid'),
    soBot: document.getElementById('so-bot'),
    soSourceInputs: document.querySelectorAll('input[name="so-mode"]'),
    soSourceListDiv: document.getElementById('so-source-list'),
    soSourceManualDiv: document.getElementById('so-source-manual'),
    soManualE: document.getElementById('so-manual-e'),
    soManualN: document.getElementById('so-manual-n'),
    soManualZ: document.getElementById('so-manual-z'),
    
    soSDist: document.getElementById('so-sdist'),
    soHr: document.getElementById('so-hr'),
    soHz: document.getElementById('so-hz'),
    soV: document.getElementById('so-v'),
    btnStakeout: document.getElementById('btn-stakeout'),
    soResult: document.getElementById('so-result'),
    soMHz: document.getElementById('so-d-hz'),
    soMDist: document.getElementById('so-d-dist'),
    soMH: document.getElementById('so-d-h'),

    // Calc
    pointListContainer: document.getElementById('point-list-container'),
    btnCalcArea: document.getElementById('btn-calc-area'),
    calcArea: document.getElementById('calc-area'),
    calcPerim: document.getElementById('calc-perim'),

    // Data
    dataTableBody: document.querySelector('#data-table tbody'),
    btnImport: document.getElementById('btn-import'),
    fileImport: document.getElementById('file-import'),
    btnExport: document.getElementById('btn-export'),
    btnClear: document.getElementById('btn-clear')
};

// --- LOGIC ---

function init() {
    console.log("App init starting...");
    loadState();
    setupEventListeners();
    renderData();
    updateStatus();
    populateSelects();
    console.log("App init done.");
}

function setupEventListeners() {
    // Navigation
    UI.tabs.forEach(tab => {
        tab.addEventListener('click', () => {
            UI.tabs.forEach(t => t.classList.remove('active'));
            UI.contents.forEach(c => c.classList.remove('active'));
            tab.classList.add('active');
            document.getElementById(tab.dataset.tab).classList.add('active');
        });
    });

    // Setup
    UI.btnSetOri.addEventListener('click', handleSetOrientation);
    
    // Setup Mode Toggles
    UI.modeKnown.addEventListener('click', () => toggleSetupMode('known'));
    UI.modeResection.addEventListener('click', () => toggleSetupMode('resection'));
    
    // Resection Events
    UI.btnAddResTarget.addEventListener('click', addResectionRow);
    UI.btnCalcRes.addEventListener('click', handleResectionCalc);
    UI.btnApplyRes.addEventListener('click', handleApplyResection);

    // Survey
    UI.btnMeasure.addEventListener('click', handleMeasurement);
    
    // Wire inputs auto-check (optional validation)
    [UI.measTop, UI.measBot].forEach(el => el.addEventListener('change', () => {
        // Could auto-calc middle or distance preview
    }));

    // Stakeout
    UI.soTargetSelect.addEventListener('change', updateStakeoutTargetInfo);
    UI.btnStakeout.addEventListener('click', handleStakeout);

    // Calc
    UI.btnCalcArea.addEventListener('click', handleAreaCalc);

    // Data
    UI.btnExport.addEventListener('click', exportCSV);
    UI.btnClear.addEventListener('click', clearData);
    
    // Global settings
    UI.kFactor.addEventListener('change', (e) => State.settings.k = parseInt(e.target.value));
    
    // Distance Mode
    UI.distModeInputs.forEach(input => {
        input.addEventListener('change', (e) => {
            if(e.target.checked) {
                State.settings.distanceMode = e.target.value;
                updateInputVisibility();
            }
        });
    });
    // Init state
    State.settings.distanceMode = 'stadia'; 
    State.settings.stakeoutSource = 'list';
    State.settings.vMode = 'zenith';

    // V-Mode Toggle
    UI.vModeInputs.forEach(input => {
        input.addEventListener('change', (e) => {
            if(e.target.checked) {
                State.settings.vMode = e.target.value;
                // Maybe update labels to say V (Z) or V (a)?
            }
        });
    });

    // Stakeout Source Toggle
    UI.soSourceInputs.forEach(input => {
        input.addEventListener('change', (e) => {
            if(e.target.checked) {
                State.settings.stakeoutSource = e.target.value;
                updateInputVisibility();
            }
        });
    });

    // Data Import
    UI.btnImport.addEventListener('click', () => UI.fileImport.click());
    UI.fileImport.addEventListener('change', handleCSVImport);
}

function updateInputVisibility() {
    const isManual = State.settings.distanceMode === 'manual';
    
    // Survey
    UI.surveyStadiaDiv.style.display = isManual ? 'none' : 'grid';
    UI.surveyManualDiv.style.display = isManual ? 'grid' : 'none';
    
    // Stakeout Readings
    UI.soStadiaDiv.style.display = isManual ? 'none' : 'grid';
    UI.soManualDiv.style.display = isManual ? 'grid' : 'none';

    // Stakeout Source
    const isSoManual = State.settings.stakeoutSource === 'manual';
    UI.soSourceListDiv.style.display = isSoManual ? 'none' : 'block';
    UI.soSourceManualDiv.style.display = isSoManual ? 'block' : 'none';
    
    // Resection - Update dynamic rows if any
    const resRows = UI.resectionList.querySelectorAll('.card');
    resRows.forEach(row => {
        const stadiaDiv = row.querySelector('.res-stadia-inputs');
        const manualDiv = row.querySelector('.res-manual-inputs');
        if(stadiaDiv && manualDiv) {
            stadiaDiv.style.display = isManual ? 'none' : 'grid';
            manualDiv.style.display = isManual ? 'grid' : 'none';
        }
    });
}

// --- HANDLERS ---

function toggleSetupMode(mode) {
    if (mode === 'known') {
        UI.setupKnown.style.display = 'block';
        UI.setupResection.style.display = 'none';
        UI.modeKnown.style.background = 'var(--primary)';
        UI.modeKnown.style.color = '#fff';
        UI.modeResection.style.background = '#ccc';
        UI.modeResection.style.color = '#333';
    } else {
        UI.setupKnown.style.display = 'none';
        UI.setupResection.style.display = 'block';
        UI.modeKnown.style.background = '#ccc';
        UI.modeKnown.style.color = '#333';
        UI.modeResection.style.background = 'var(--primary)';
        UI.modeResection.style.color = '#fff';
        
        // Populate Selects in rows if needed, or if first time add a row
        if(UI.resectionList.children.length === 0) {
            addResectionRow();
            addResectionRow(); // Start with 2
        }
    }
}

function addResectionRow() {
    const div = document.createElement('div');
    div.className = 'card';
    div.style.padding = '10px';
    div.style.background = '#f9f9f9';
    
    // Build options via DOM for safety
    const select = document.createElement('select');
    select.className = 'res-target';
    
    const defOpt = document.createElement('option');
    defOpt.value = "";
    defOpt.textContent = "Select Target...";
    select.appendChild(defOpt);
    
    State.points.forEach(p => {
        const opt = document.createElement('option');
        opt.value = p.id;
        opt.textContent = p.id;
        select.appendChild(opt);
    });

    div.innerHTML = `
        <div class="form-group">
            <label>Target Point</label>
            <div class="select-container"></div> 
        </div>
        
        <div class="input-grid three-col res-stadia-inputs" style="display: ${State.settings.distanceMode === 'manual' ? 'none' : 'grid'}">
            <div class="form-group"><label>Top</label><input type="number" class="res-top" step="0.001"></div>
            <div class="form-group"><label>Mid</label><input type="number" class="res-mid" step="0.001"></div>
            <div class="form-group"><label>Bot</label><input type="number" class="res-bot" step="0.001"></div>
        </div>

        <div class="input-grid res-manual-inputs" style="display: ${State.settings.distanceMode === 'manual' ? 'grid' : 'none'}">
            <div class="form-group"><label>Slope Dist</label><input type="number" class="res-sdist" step="0.001"></div>
            <div class="form-group"><label>Tgt Ht</label><input type="number" class="res-hr" step="0.001"></div>
        </div>

        <div class="input-grid">
            <div class="form-group">
                <label>Hz (g)</label>
                <input type="number" class="res-hz" step="0.0001">
            </div>
            <div class="form-group">
                <label>V (g)</label>
                <input type="number" class="res-v" step="0.0001" placeholder="100.00">
            </div>
        </div>
        <button class="action-btn sm danger remove-row" style="margin-top:5px;">Remove</button>
    `;
    
    // Inject select
    div.querySelector('.select-container').appendChild(select);
    
    div.querySelector('.remove-row').addEventListener('click', () => {
        div.remove();
    });
    
    UI.resectionList.appendChild(div);
}

function handleResectionCalc() {
    const rows = UI.resectionList.querySelectorAll('.card');
    if (rows.length < 2) {
        alert("Need at least 2 measurements for resection.");
        return;
    }

    const obs = [];
    const k = State.settings.k;
    let error = false;

    rows.forEach(row => {
        const pid = row.querySelector('.res-target').value;
        const hz = parseFloat(row.querySelector('.res-hz').value);
        const vInput = row.querySelector('.res-v').value;
        const v = vInput ? parseFloat(vInput) : 100.0;
        
        if (!pid || isNaN(hz)) return;

        const pt = State.points.find(p => p.id === pid);
        if(!pt) return;

        let horizDist = 0;

        if (State.settings.distanceMode === 'manual') {
            // Manual Mode: Slope Dist + Target Height
            const sDist = parseFloat(row.querySelector('.res-sdist').value);
            // We don't need Tgt Height for Resection 2D pos, only for Z. 
            // But if we wanted to calculate Reduced distance from Slope, we need V.
            // Horizontal Distance = Slope * sin(Z_rad)
            if (isNaN(sDist)) return;
            
            const zRad = MathUtils.gradsToRad(v);
            horizDist = sDist * Math.sin(zRad);
            
        } else {
            // Stadia Mode
            const top = parseFloat(row.querySelector('.res-top').value);
            const mid = parseFloat(row.querySelector('.res-mid').value);
            const bot = parseFloat(row.querySelector('.res-bot').value);
            if (isNaN(top) || isNaN(bot)) return;
            
            horizDist = MathUtils.calcDist(top, bot, v, k);
        }
        
        obs.push({
            e: pt.e,
            n: pt.n,
            dist: horizDist,
            hz: hz
        });
    });

    if (obs.length < 2) {
        alert("Please enter at least 2 valid observations.");
        return;
    }

    const result = MathUtils.solveResection(obs);
    
    if (result) {
        UI.resStnE.textContent = result.e.toFixed(3);
        UI.resStnN.textContent = result.n.toFixed(3);
        UI.resOri.textContent = result.ori.toFixed(4);
        UI.resResult.style.display = 'block';
        
        // Store temp result
        State.tempResection = result;
    } else {
        alert("Resection calculation failed (Singular matrix or no convergence). Check inputs.");
    }
}

function handleApplyResection() {
    if (!State.tempResection) return;
    
    const { e, n, ori } = State.tempResection;
    const hi = parseFloat(UI.resHi.value) || 0;
    
    // Note: Resection usually solves 2D. Height is separate.
    // If we want Z, we need Z of targets + Delta H.
    // Z_stn = Z_tgt - dH + hi - mid
    // Simple average of Z from all targets?
    // Let's keep Z as 0 or user input? 
    // The current UI shows E, N. I'll leave Z as 0 (or previous value).
    
    State.station = { e, n, z: 0, hi, set: true };
    State.orientation = { azimuthCorrection: ori, set: true };
    
    saveState();
    updateStatus();
    
    // Update inputs to reflect
    UI.stnE.value = e.toFixed(3);
    UI.stnN.value = n.toFixed(3);
    
    alert("Station and Orientation updated from Resection.");
}

function handleSetOrientation() {
    console.log("Setting Orientation...");
    const e = parseFloat(UI.stnE.value);
    const n = parseFloat(UI.stnN.value);
    const z = parseFloat(UI.stnZ.value) || 0;
    const hi = parseFloat(UI.stnHi.value) || 0;
    
    const bsE = parseFloat(UI.bsE.value);
    const bsN = parseFloat(UI.bsN.value);
    const bsHz = parseFloat(UI.bsHz.value) || 0;

    if (isNaN(e) || isNaN(n) || isNaN(bsE) || isNaN(bsN)) {
        alert("Please enter valid coordinates for Station and Backsight.");
        return;
    }

    // 1. Calc Azimuth Station -> Backsight
    const trueAzimuth = MathUtils.azimuth(e, n, bsE, bsN);

    // 2. Calc Orientation Correction (Angle to add to readings to get Azimuth)
    // Azimuth = Reading + Correction => Correction = Azimuth - Reading
    let correction = trueAzimuth - bsHz;
    correction = MathUtils.normalizeGrads(correction);

    State.station = { e, n, z, hi, set: true };
    State.orientation = { azimuthCorrection: correction, set: true };

    saveState();
    updateStatus();
    alert(`Station Set!\nAzimuth to BS: ${trueAzimuth.toFixed(4)}g\nOrientation Cor: ${correction.toFixed(4)}g`);
}

function handleMeasurement() {
    if (!State.station.set || !State.orientation.set) {
        alert("Please set up Station and Orientation first.");
        return;
    }

    const pid = UI.measPid.value || ("PT" + (State.points.length + 1));
    const hz = parseFloat(UI.measHz.value);
    const v = parseFloat(UI.measV.value);

    if (isNaN(hz) || isNaN(v)) {
        alert("Please fill angles.");
        return;
    }

    let horizDist, dH_component, targetHeight;
    const k = State.settings.k;

    if (State.settings.distanceMode === 'manual') {
        const sDist = parseFloat(UI.measSDist.value);
        const hr = parseFloat(UI.measHr.value);
        if (isNaN(sDist) || isNaN(hr)) {
            alert("Please fill Slope Distance and Target Height.");
            return;
        }
        
        // Manual Calc
        let z = v;
        if(State.settings.vMode === 'alpha') {
            z = 100 - v;
        }
        
        const zRad = MathUtils.gradsToRad(z);
        horizDist = sDist * Math.sin(zRad);
        // Vertical Diff (Inst Axis to Target Axis) = sDist * cos(Z)
        // dH_component (Inst Axis to Target Center)
        dH_component = sDist * Math.cos(zRad);
        targetHeight = hr;

    } else {
        // Stadia Calc
        const top = parseFloat(UI.measTop.value);
        const mid = parseFloat(UI.measMid.value);
        const bot = parseFloat(UI.measBot.value);

        if (isNaN(top) || isNaN(bot) || isNaN(mid)) {
            alert("Please fill wire readings.");
            return;
        }

        horizDist = MathUtils.calcDist(top, bot, v, k);
        dH_component = MathUtils.calcVertDiffComponent(top, bot, v, k);
        targetHeight = mid;
    }

    // 2. Calc Coordinates
    // Azimuth = Measured Hz + Correction
    const azimuth = MathUtils.normalizeGrads(hz + State.orientation.azimuthCorrection);
    const coord2D = MathUtils.polar(State.station.e, State.station.n, azimuth, horizDist);

    // 3. Calc Height
    // Target Z = Station Z + InstHeight + dH_component - TargetHeight
    const targetZ = State.station.z + State.station.hi + dH_component - targetHeight;

    const point = {
        id: pid,
        e: coord2D.e,
        n: coord2D.n,
        z: targetZ,
        code: 'MEAS',
        ts: new Date().toISOString()
    };

    // Save
    State.points.push(point);
    saveState();
    renderData();
    populateSelects();

    // Show Result
    UI.resDist.textContent = horizDist.toFixed(3);
    UI.resE.textContent = point.e.toFixed(3);
    UI.resN.textContent = point.n.toFixed(3);
    UI.resH.textContent = point.z.toFixed(3);
    UI.resCard.style.display = 'block';

    // Auto increment ID logic could go here
}

function handleStakeout() {
    console.log("handleStakeout called");
    if (!State.station.set || !State.orientation.set) {
        console.log("Station not set");
        alert("Please set up Station and Orientation first.");
        return;
    }

    let target;

    if (State.settings.stakeoutSource === 'manual') {
        const e = parseFloat(UI.soManualE.value);
        const n = parseFloat(UI.soManualN.value);
        const z = parseFloat(UI.soManualZ.value);
        if (isNaN(e) || isNaN(n) || isNaN(z)) {
            alert("Please fill Target Coordinates.");
            return;
        }
        target = { e, n, z };
    } else {
        const targetId = UI.soTargetSelect.value;
        if (!targetId) return;
        target = State.points.find(p => p.id === targetId);
        if (!target) return;
    }

    const hz = parseFloat(UI.soHz.value);
    const v = parseFloat(UI.soV.value);

    if (isNaN(hz) || isNaN(v)) {
        alert("Please fill angles.");
        return;
    }

    let horizDistCurrent, dH_component, targetHeight;
    const k = State.settings.k;

    if (State.settings.distanceMode === 'manual') {
        const sDist = parseFloat(UI.soSDist.value);
        const hr = parseFloat(UI.soHr.value);
        if (isNaN(sDist) || isNaN(hr)) {
            alert("Please fill Slope Distance and Target Height.");
            return;
        }
        let z = v;
        if(State.settings.vMode === 'alpha') {
            z = 100 - v;
        }
        
        const zRad = MathUtils.gradsToRad(z);
        horizDistCurrent = sDist * Math.sin(zRad);
        dH_component = sDist * Math.cos(zRad);
        targetHeight = hr;
    } else {
        const top = parseFloat(UI.soTop.value);
        const mid = parseFloat(UI.soMid.value);
        const bot = parseFloat(UI.soBot.value);
        if (isNaN(top) || isNaN(bot) || isNaN(mid)) {
            alert("Please fill wire readings.");
            return;
        }

        horizDistCurrent = MathUtils.calcDist(top, bot, v, k);
        dH_component = MathUtils.calcVertDiffComponent(top, bot, v, k);
        targetHeight = mid;
    }

    // Current position calculation
    const azimuthCurrent = MathUtils.normalizeGrads(hz + State.orientation.azimuthCorrection);
    const zCurrent = State.station.z + State.station.hi + dH_component - targetHeight;
    // const currentPos = MathUtils.polar(State.station.e, State.station.n, azimuthCurrent, horizDistCurrent); 
    // ^ Don't really need full coords for stakeout deltas unless we want to display them

    // Required calculations
    const reqAzimuth = MathUtils.azimuth(State.station.e, State.station.n, target.e, target.n);
    const reqDist = MathUtils.dist2D(State.station.e, State.station.n, target.e, target.n);
    
    // Deltas
    // Turn: Right + / Left -
    let dHz = reqAzimuth - azimuthCurrent; // If positive, turn right (increase angle). If negative, turn left.
    // Normalize to -200 to +200 for easier reading
    if (dHz > 200) dHz -= 400;
    if (dHz < -200) dHz += 400;
    
    const dDist = reqDist - horizDistCurrent; // Positive = Move Back (Away), Negative = Move Forward (Closer)
    const dH = target.z - zCurrent; // Positive = Fill (Go up), Negative = Cut (Go down)

    // Display
    UI.soMHz.textContent = `${dHz >= 0 ? 'Right' : 'Left'} ${Math.abs(dHz).toFixed(4)}g`;
    UI.soMDist.textContent = `${dDist >= 0 ? 'Back' : 'Fwd'} ${Math.abs(dDist).toFixed(3)}m`;
    UI.soMH.textContent = `${dH >= 0 ? 'Fill' : 'Cut'} ${Math.abs(dH).toFixed(3)}m`;
    
    UI.soResult.style.display = 'block';
}

function handleAreaCalc() {
    const checkboxes = document.querySelectorAll('.pt-select:checked');
    const selectedPoints = Array.from(checkboxes).map(cb => {
        return State.points.find(p => p.id === cb.value);
    });

    if (selectedPoints.length < 3) {
        alert("Select at least 3 points.");
        return;
    }

    const area = MathUtils.calcPolygonArea(selectedPoints);
    const perim = MathUtils.calcPerimeter(selectedPoints);

    UI.calcArea.textContent = area.toFixed(3);
    UI.calcPerim.textContent = perim.toFixed(3);
}

// --- DATA MANAGEMENT ---

function saveState() {
    localStorage.setItem('webts_state', JSON.stringify(State));
}

function loadState() {
    const saved = localStorage.getItem('webts_state');
    if (saved) {
        const parsed = JSON.parse(saved);
        State.station = parsed.station || State.station;
        State.orientation = parsed.orientation || State.orientation;
        State.points = parsed.points || [];
        State.settings = parsed.settings || { k: 100 };
        
        // Restore inputs
        UI.kFactor.value = State.settings.k;
        if (State.station.set) {
            UI.stnE.value = State.station.e;
            UI.stnN.value = State.station.n;
            UI.stnZ.value = State.station.z;
            UI.stnHi.value = State.station.hi;
        }
    }
}

function updateStatus() {
    UI.statusStn.textContent = State.station.set ? "STN: Set" : "STN: Not Set";
    UI.statusStn.style.color = State.station.set ? "#4caf50" : "#aaa";
    
    UI.statusOri.textContent = State.orientation.set ? "ORI: Set" : "ORI: Not Set";
    UI.statusOri.style.color = State.orientation.set ? "#4caf50" : "#aaa";
}

function renderData() {
    UI.dataTableBody.innerHTML = '';
    
    // Also update calc list
    UI.pointListContainer.innerHTML = '';

    State.points.forEach(p => {
        // Table Row - Use textContent for safety
        const tr = document.createElement('tr');
        
        const tdId = document.createElement('td'); tdId.textContent = p.id;
        const tdE = document.createElement('td'); tdE.textContent = p.e.toFixed(3);
        const tdN = document.createElement('td'); tdN.textContent = p.n.toFixed(3);
        const tdZ = document.createElement('td'); tdZ.textContent = p.z.toFixed(3);
        const tdC = document.createElement('td'); tdC.textContent = p.code;
        
        tr.append(tdId, tdE, tdN, tdZ, tdC);
        UI.dataTableBody.appendChild(tr);

        // Calc List Item
        const label = document.createElement('label');
        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.className = 'pt-select';
        checkbox.value = p.id; // Value attribute is safe
        
        label.appendChild(checkbox);
        label.appendChild(document.createTextNode(" " + p.id));
        
        UI.pointListContainer.appendChild(label);
    });
}

function populateSelects() {
    // Stakeout target
    UI.soTargetSelect.innerHTML = '<option value="">Select Point...</option>';
    State.points.forEach(p => {
        const opt = document.createElement('option');
        opt.value = p.id;
        opt.textContent = p.id;
        UI.soTargetSelect.appendChild(opt);
    });
}

function updateStakeoutTargetInfo() {
    const pid = UI.soTargetSelect.value;
    const p = State.points.find(x => x.id === pid);
    if (p) {
        UI.soTgtE.textContent = p.e.toFixed(3);
        UI.soTgtN.textContent = p.n.toFixed(3);
        UI.soTargetInfo.style.display = 'block';
    } else {
        UI.soTargetInfo.style.display = 'none';
    }
}

function clearData() {
    if(confirm("Are you sure you want to clear all data? This cannot be undone.")) {
        State.points = [];
        State.station.set = false;
        State.orientation.set = false;
        saveState();
        renderData();
        updateStatus();
        populateSelects();
    }
}

function exportCSV() {
    let csv = "PointID,East,North,Elevation,Code\n";
    State.points.forEach(p => {
        csv += `${p.id},${p.e},${p.n},${p.z},${p.code}\n`;
    });
    
    const blob = new Blob([csv], { type: 'text/csv' });
    const url = window.URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.setAttribute('hidden', '');
    a.setAttribute('href', url);
    a.setAttribute('download', 'survey_data.csv');
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
}

function handleCSVImport(event) {
    const file = event.target.files[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
        const text = e.target.result;
        const lines = text.split('\n');
        
        let importedCount = 0;
        
        // Determine format based on header
        const header = lines[0].toLowerCase();
        // Support: Comma or Tab
        const sep = header.includes('\t') ? '\t' : ',';
        
        let map = { id: 0, e: 1, n: 2, z: 3, c: 4 }; // Default standard
        const startIndex = 1; // Always skip header for safety if we detect one

        if (header.includes('id') || header.includes('point')) {
            const hCols = header.split(sep).map(s => s.trim());
            
            const idxId = hCols.findIndex(c => c === 'id' || c === 'pointid' || c === 'point');
            const idxE = hCols.findIndex(c => c === 'east' || c === 'e' || c === 'y'); // Y is often East in CAD/Math but user said Y=Easting
            const idxN = hCols.findIndex(c => c === 'north' || c === 'n' || c === 'x'); // X is often North in Surveying (local grids sometimes)
            const idxZ = hCols.findIndex(c => c === 'elev' || c === 'z' || c === 'h' || c === 'height');
            const idxC = hCols.findIndex(c => c === 'code' || c === 'cd');

            if (idxId >= 0 && idxE >= 0 && idxN >= 0) {
                map = { id: idxId, e: idxE, n: idxN, z: idxZ, c: idxC };
            }
        }

        for (let i = startIndex; i < lines.length; i++) {
            const line = lines[i].trim();
            if (!line) continue;
            
            const cols = line.split(sep);
            if (cols.length < 2) continue; 
            
            const id = cols[map.id].trim();
            const eVal = parseFloat(cols[map.e]);
            const nVal = parseFloat(cols[map.n]);
            const zVal = (map.z >= 0) ? (parseFloat(cols[map.z]) || 0) : 0;
            const code = (map.c >= 0 && cols[map.c]) ? cols[map.c].trim() : '';
            
            if (isNaN(eVal) || isNaN(nVal)) continue;
            
            // Check if exists
            const existingIndex = State.points.findIndex(p => p.id === id);
            const point = { id, e: eVal, n: nVal, z: zVal, code, ts: new Date().toISOString() };
            
            if (existingIndex >= 0) {
                State.points[existingIndex] = point; // Overwrite
            } else {
                State.points.push(point);
            }
            importedCount++;
        }
        
        saveState();
        renderData();
        populateSelects();
        alert(`Imported ${importedCount} points.`);
        UI.fileImport.value = ''; // Reset input
    };
    reader.readAsText(file);
}

// Start
init();
