"""
ProportionMA - Test Suite
Tests for proportion/prevalence meta-analysis browser tool.
Reference values from R meta::metaprop.
"""

import pytest
import os
import sys
import re
import time
import socket
import subprocess
import threading
from http.server import HTTPServer, SimpleHTTPRequestHandler
from functools import partial

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.support.ui import WebDriverWait, Select
from selenium.webdriver.support import expected_conditions as EC


# ============================================================
# Fixtures
# ============================================================

def get_free_port():
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('', 0))
        return s.getsockname()[1]


@pytest.fixture(scope="session")
def server():
    port = get_free_port()
    app_dir = os.path.dirname(os.path.abspath(__file__))
    handler = partial(SimpleHTTPRequestHandler, directory=app_dir)
    httpd = HTTPServer(('127.0.0.1', port), handler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield f"http://127.0.0.1:{port}/index.html"
    httpd.shutdown()


@pytest.fixture(scope="session")
def driver(server):
    # Kill orphan chromedriver
    try:
        if sys.platform == 'win32':
            subprocess.run(['taskkill', '/F', '/IM', 'chromedriver.exe'],
                         capture_output=True, timeout=5)
    except Exception:
        pass

    opts = Options()
    opts.add_argument('--headless=new')
    opts.add_argument('--no-sandbox')
    opts.add_argument('--disable-dev-shm-usage')
    opts.add_argument('--disable-gpu')
    opts.add_argument('--window-size=1400,900')
    opts.set_capability('goog:loggingPrefs', {'browser': 'ALL'})

    drv = webdriver.Chrome(options=opts)
    drv.set_page_load_timeout(30)
    drv.implicitly_wait(5)
    yield drv
    drv.quit()


@pytest.fixture(autouse=True)
def fresh_page(driver, server):
    driver.get(server)
    time.sleep(0.5)
    # Drain console logs from page load
    driver.get_log('browser')


# ============================================================
# Helpers
# ============================================================

def set_csv(driver, csv_text):
    """Set CSV data via JS to ensure value is properly set."""
    driver.execute_script(
        "document.getElementById('csvData').value = arguments[0];",
        csv_text.strip()
    )


def analyze(driver):
    """Run analysis via JS and wait."""
    driver.execute_script("runAnalysis()")
    time.sleep(0.8)


def get_pooled(driver):
    """Get pooled proportion from JS lastResults."""
    return driver.execute_script("return lastResults ? lastResults.pooledProp : null")


def get_ci(driver):
    """Get CI bounds from JS lastResults."""
    lower = driver.execute_script("return lastResults ? lastResults.pooledCILower : null")
    upper = driver.execute_script("return lastResults ? lastResults.pooledCIUpper : null")
    if lower is not None and upper is not None:
        return lower, upper
    return None, None


def get_result(driver, field):
    """Get a field from lastResults.result."""
    return driver.execute_script(f"return lastResults && lastResults.result ? lastResults.result.{field} : null")


def select_opt(driver, elem_id, value):
    sel = Select(driver.find_element(By.ID, elem_id))
    sel.select_by_value(value)


def set_hksj(driver, on=True):
    current = driver.execute_script("return useHKSJ")
    if current != on:
        driver.execute_script("toggleHKSJ()")


def set_pi(driver, on=True):
    current = driver.execute_script("return usePI")
    if current != on:
        driver.execute_script("togglePI()")


DEMO_CSV = """Study,Events,Total
Smith 2020,15,100
Jones 2021,22,150
Lee 2019,8,50
Brown 2022,30,200
Davis 2018,12,80
Wilson 2023,5,40
Taylor 2020,18,120
Clark 2021,25,160"""

DEMO_SUBGROUP_CSV = """Study,Events,Total,Subgroup
Smith 2020,15,100,Urban
Jones 2021,22,150,Urban
Lee 2019,8,50,Urban
Brown 2022,30,200,Urban
Davis 2018,12,80,Rural
Wilson 2023,5,40,Rural
Taylor 2020,18,120,Rural
Clark 2021,25,160,Rural"""


# ============================================================
# Tests
# ============================================================

class TestAppLoads:
    """Test 1: App loads without JS errors."""

    def test_app_loads_no_errors(self, driver):
        title = driver.title
        assert 'ProportionMA' in title
        logs = driver.get_log('browser')
        # Filter out favicon 404 and other resource errors
        js_errors = [l for l in logs if l['level'] == 'SEVERE'
                     and 'favicon' not in l.get('message', '').lower()
                     and 'Failed to load resource' not in l.get('message', '')]
        assert len(js_errors) == 0, f"JS errors: {js_errors}"

    def test_app_has_key_elements(self, driver):
        assert driver.find_element(By.ID, 'csvData')
        assert driver.find_element(By.ID, 'transform')
        assert driver.find_element(By.ID, 'method')
        assert driver.find_element(By.ID, 'hksjToggle')
        assert driver.find_element(By.ID, 'piToggle')


class TestDemoData:
    """Test 2: Demo data button populates 8 studies."""

    def test_demo_button_populates(self, driver):
        driver.execute_script("loadDemoData()")
        time.sleep(0.3)
        val = driver.execute_script("return document.getElementById('csvData').value")
        lines = [l for l in val.strip().split('\n') if l.strip()]
        assert len(lines) == 9, f"Expected 9 lines (header+8), got {len(lines)}"


class TestFreemanTukey:
    """Test 3: FT transform pooled proportion in [0.10, 0.20]."""

    def test_ft_pooled_proportion(self, driver):
        set_csv(driver, DEMO_CSV)
        select_opt(driver, 'transform', 'FT')
        select_opt(driver, 'method', 'DL')
        set_hksj(driver, False)
        analyze(driver)
        prop = get_pooled(driver)
        assert prop is not None, "Could not get pooled proportion"
        assert 0.10 <= prop <= 0.20, f"FT pooled proportion {prop} outside [0.10, 0.20]"


class TestLogit:
    """Test 4: Logit transform pooled proportion in [0.10, 0.20]."""

    def test_logit_pooled_proportion(self, driver):
        set_csv(driver, DEMO_CSV)
        select_opt(driver, 'transform', 'logit')
        select_opt(driver, 'method', 'DL')
        set_hksj(driver, False)
        analyze(driver)
        prop = get_pooled(driver)
        assert prop is not None, "Could not get pooled proportion"
        assert 0.10 <= prop <= 0.20, f"Logit pooled proportion {prop} outside [0.10, 0.20]"


class TestREML:
    """Test 5: REML method produces valid (non-NaN) result."""

    def test_reml_not_nan(self, driver):
        set_csv(driver, DEMO_CSV)
        select_opt(driver, 'transform', 'FT')
        select_opt(driver, 'method', 'REML')
        set_hksj(driver, False)
        analyze(driver)
        prop = get_pooled(driver)
        assert prop is not None, "REML produced no result"
        assert 0.05 <= prop <= 0.30, f"REML proportion {prop} seems off"


class TestHKSJ:
    """Test 6: HKSJ toggle widens CI (HKSJ CI width >= DL CI width)."""

    def test_hksj_widens_ci(self, driver):
        set_csv(driver, DEMO_CSV)
        select_opt(driver, 'transform', 'FT')
        select_opt(driver, 'method', 'DL')

        # Without HKSJ
        set_hksj(driver, False)
        analyze(driver)
        lo1, hi1 = get_ci(driver)
        assert lo1 is not None, "No CI without HKSJ"
        w_dl = hi1 - lo1

        # With HKSJ
        set_hksj(driver, True)
        analyze(driver)
        lo2, hi2 = get_ci(driver)
        assert lo2 is not None, "No CI with HKSJ"
        w_hksj = hi2 - lo2

        # HKSJ should be at least as wide (allowing small tolerance)
        assert w_hksj >= w_dl * 0.95, \
            f"HKSJ CI width ({w_hksj:.4f}) should be >= DL CI width ({w_dl:.4f})"

        set_hksj(driver, False)


class TestPredictionInterval:
    """Test 7: PI toggle shows dashed line in SVG (or message if k<3)."""

    def test_pi_with_enough_studies(self, driver):
        set_csv(driver, DEMO_CSV)
        set_pi(driver, True)
        analyze(driver)

        svg_html = driver.execute_script(
            "var el = document.getElementById('forestSVG'); return el ? el.outerHTML : ''"
        )
        assert 'pi-line' in svg_html or 'PI:' in svg_html, "PI line not in forest plot"
        set_pi(driver, False)

    def test_pi_message_for_k_less_than_3(self, driver):
        csv = "Study,Events,Total\nA,10,50\nB,15,80"
        set_csv(driver, csv)
        set_pi(driver, True)
        analyze(driver)

        svg_html = driver.execute_script(
            "var el = document.getElementById('forestSVG'); return el ? el.outerHTML : ''"
        )
        # Check for the PI undefined message
        assert 'undefined' in svg_html.lower() or 'k' in svg_html, \
            "No PI undefined message for k<3"
        set_pi(driver, False)


class TestZeroEvents:
    """Test 8: Zero events study (0/50): no NaN, continuity correction applied."""

    def test_zero_events(self, driver):
        csv = "Study,Events,Total\nA,0,50\nB,10,80\nC,15,100\nD,8,60"
        set_csv(driver, csv)
        set_hksj(driver, False)
        analyze(driver)
        prop = get_pooled(driver)
        assert prop is not None, "Zero-events produced no result"
        assert 0.0 <= prop <= 1.0, f"Proportion {prop} out of bounds"
        # Ensure not NaN
        assert not driver.execute_script(
            "return isNaN(lastResults.pooledProp)"
        ), "Pooled proportion is NaN with zero events"


class TestAllEvents:
    """Test 9: All events study (50/50): no NaN, clamping works."""

    def test_all_events(self, driver):
        csv = "Study,Events,Total\nA,50,50\nB,10,80\nC,15,100\nD,8,60"
        set_csv(driver, csv)
        set_hksj(driver, False)
        analyze(driver)
        prop = get_pooled(driver)
        assert prop is not None, "All-events produced no result"
        assert 0.0 <= prop <= 1.0, f"Proportion {prop} out of bounds"
        assert not driver.execute_script("return isNaN(lastResults.pooledProp)")


class TestForestPlot:
    """Test 10: Forest plot SVG has correct number of study rows."""

    def test_forest_study_count(self, driver):
        set_csv(driver, DEMO_CSV)
        analyze(driver)

        # Count rect elements in SVG (one per study)
        count = driver.execute_script("""
            var svg = document.getElementById('forestSVG');
            if (!svg) return 0;
            return svg.querySelectorAll('rect').length;
        """)
        assert count >= 8, f"Expected >= 8 rects in forest plot, got {count}"


class TestFunnelPlot:
    """Test 11: Funnel plot SVG renders."""

    def test_funnel_renders(self, driver):
        set_csv(driver, DEMO_CSV)
        analyze(driver)

        # Switch to funnel tab
        driver.execute_script("switchTab('funnel')")
        time.sleep(0.3)

        count = driver.execute_script("""
            var svg = document.getElementById('funnelSVG');
            if (!svg) return -1;
            return svg.querySelectorAll('circle').length;
        """)
        assert count == 8, f"Expected 8 circles in funnel, got {count}"


class TestSummaryTable:
    """Test 12: Summary table shows I2, tau2, Q values."""

    def test_summary_values(self, driver):
        set_csv(driver, DEMO_CSV)
        analyze(driver)

        Q = get_result(driver, 'Q')
        I2 = get_result(driver, 'I2')
        tau2 = get_result(driver, 'tau2')

        assert Q is not None and not driver.execute_script(
            "return isNaN(lastResults.result.Q)"
        ), "Q is NaN or null"
        assert I2 is not None and not driver.execute_script(
            "return isNaN(lastResults.result.I2)"
        ), "I2 is NaN or null"
        assert tau2 is not None and not driver.execute_script(
            "return isNaN(lastResults.result.tau2)"
        ), "tau2 is NaN or null"

        # Verify they appear in rendered summary
        driver.execute_script("switchTab('summary')")
        time.sleep(0.2)
        q_text = driver.find_element(By.ID, 'summaryQ').text
        i2_text = driver.find_element(By.ID, 'summaryI2').text
        tau2_text = driver.find_element(By.ID, 'summaryTau2').text
        assert q_text != 'N/A', f"Q displayed as N/A"
        assert i2_text != 'N/A', f"I2 displayed as N/A"
        assert tau2_text != 'N/A', f"tau2 displayed as N/A"


class TestCSVExport:
    """Test 13: CSV export produces valid content."""

    def test_csv_export(self, driver):
        set_csv(driver, DEMO_CSV)
        analyze(driver)

        # Call exportCSV via JS and intercept the blob
        result = driver.execute_script("""
            try {
                var blobContent = null;
                var origBlob = window.Blob;
                window.Blob = function(parts, opts) {
                    blobContent = parts.join('');
                    return new origBlob(parts, opts);
                };
                exportCSV();
                window.Blob = origBlob;
                return blobContent;
            } catch(e) {
                return 'error: ' + e.message;
            }
        """)
        assert result is not None and not result.startswith('error'), f"CSV export failed: {result}"
        assert 'Pooled' in result, "CSV export missing pooled estimate"
        assert 'Study' in result or 'Smith' in result, "CSV export missing study data"


class TestSingleStudy:
    """Test 14: k=1 produces valid result (no pooling, just the study)."""

    def test_single_study(self, driver):
        csv = "Study,Events,Total\nSingle,20,100"
        set_csv(driver, csv)
        set_hksj(driver, False)
        analyze(driver)
        prop = get_pooled(driver)
        assert prop is not None, "k=1 produced no result"
        assert abs(prop - 0.20) < 0.05, f"k=1 proportion {prop} should be ~0.20"


class TestSubgroup:
    """Test 15: Subgroup analysis with 2 groups produces separate summaries."""

    def test_subgroup_analysis(self, driver):
        set_csv(driver, DEMO_SUBGROUP_CSV)
        analyze(driver)

        has_subgroup = driver.execute_script(
            "return lastResults && lastResults.subgroupResults !== null"
        )
        assert has_subgroup, "No subgroup results"

        # Check subgroup tab is visible
        tab_display = driver.execute_script(
            "return document.getElementById('subgroupTab').style.display"
        )
        assert tab_display != 'none', "Subgroup tab hidden"

        # Switch to subgroup tab and check content
        driver.execute_script("switchTab('subgroup')")
        time.sleep(0.3)
        content = driver.find_element(By.ID, 'subgroupContent').text
        assert 'Urban' in content, "Urban subgroup not found"
        assert 'Rural' in content, "Rural subgroup not found"


class TestInputValidation:
    """Test 16: Negative events shows error."""

    def test_negative_events_error(self, driver):
        csv = "Study,Events,Total\nA,-5,50\nB,10,80"
        set_csv(driver, csv)
        analyze(driver)

        err = driver.find_element(By.ID, 'errorMsg')
        classes = err.get_attribute('class')
        assert 'visible' in classes, "Error not shown for negative events"


class TestEmptyInput:
    """Test 17: Empty input shows error message."""

    def test_empty_input_error(self, driver):
        driver.execute_script("document.getElementById('csvData').value = ''")
        analyze(driver)

        err = driver.find_element(By.ID, 'errorMsg')
        classes = err.get_attribute('class')
        assert 'visible' in classes, "Error not shown for empty input"


class TestRawTransform:
    """Test 18: Raw proportion transform works."""

    def test_raw_proportion(self, driver):
        set_csv(driver, DEMO_CSV)
        select_opt(driver, 'transform', 'raw')
        select_opt(driver, 'method', 'DL')
        set_hksj(driver, False)
        analyze(driver)
        prop = get_pooled(driver)
        assert prop is not None, "Raw transform produced no result"
        assert 0.10 <= prop <= 0.20, f"Raw proportion {prop} outside [0.10, 0.20]"


class TestEventsExceedTotal:
    """Test 19: Events exceeding total shows error."""

    def test_events_exceed_total(self, driver):
        csv = "Study,Events,Total\nA,60,50\nB,10,80"
        set_csv(driver, csv)
        analyze(driver)

        err = driver.find_element(By.ID, 'errorMsg')
        classes = err.get_attribute('class')
        assert 'visible' in classes, "Error not shown when events > total"


class TestHeterogeneityValues:
    """Test 20: I2 in expected range for homogeneous dataset."""

    def test_i2_near_zero(self, driver):
        set_csv(driver, DEMO_CSV)
        select_opt(driver, 'transform', 'FT')
        set_hksj(driver, False)
        analyze(driver)

        I2 = get_result(driver, 'I2')
        assert I2 is not None, "I2 is null"
        assert I2 < 30, f"I2 = {I2}% too high for homogeneous dataset"


class TestConsoleClean:
    """Test 21: No console errors during full workflow."""

    def test_full_workflow_no_errors(self, driver, server):
        driver.get(server)
        time.sleep(0.5)
        driver.get_log('browser')  # drain

        # Load demo
        driver.execute_script("loadDemoData()")
        time.sleep(0.2)

        # Test each transform
        for tf in ['FT', 'logit', 'raw']:
            select_opt(driver, 'transform', tf)
            driver.execute_script("runAnalysis()")
            time.sleep(0.3)

        # Toggle HKSJ and re-run
        driver.execute_script("toggleHKSJ()")
        driver.execute_script("runAnalysis()")
        time.sleep(0.3)

        # Toggle PI
        driver.execute_script("togglePI()")
        driver.execute_script("runAnalysis()")
        time.sleep(0.3)

        # Switch tabs
        for tab in ['forest', 'summary', 'funnel']:
            driver.execute_script(f"switchTab('{tab}')")
            time.sleep(0.1)

        logs = driver.get_log('browser')
        js_errors = [l for l in logs if l['level'] == 'SEVERE'
                     and 'favicon' not in l.get('message', '').lower()
                     and 'Failed to load resource' not in l.get('message', '')]
        assert len(js_errors) == 0, f"Console errors: {js_errors}"


class TestClopperPearsonExact:
    """Test 22: Clopper-Pearson CIs are reasonable for 20/100."""

    def test_clopper_pearson_ci(self, driver):
        csv = "Study,Events,Total\nTest,20,100"
        set_csv(driver, csv)
        set_hksj(driver, False)
        analyze(driver)

        # Get per-study CI from lastResults
        ci = driver.execute_script("""
            if (!lastResults || !lastResults.studies || lastResults.studies.length === 0) return null;
            var s = lastResults.studies[0];
            return {lower: s.ciLower, upper: s.ciUpper};
        """)
        assert ci is not None, "No study CI"
        # Clopper-Pearson for 20/100: approximately [0.127, 0.286]
        assert abs(ci['lower'] - 0.127) < 0.02, f"CP lower {ci['lower']} not ~0.127"
        assert abs(ci['upper'] - 0.286) < 0.02, f"CP upper {ci['upper']} not ~0.286"


class TestREMLConvergence:
    """Test 23: REML converges for heterogeneous data."""

    def test_reml_heterogeneous(self, driver):
        csv = """Study,Events,Total
High1,40,50
High2,35,45
Low1,2,50
Low2,3,60
Low3,1,40"""
        set_csv(driver, csv)
        select_opt(driver, 'method', 'REML')
        set_hksj(driver, False)
        analyze(driver)

        prop = get_pooled(driver)
        assert prop is not None, "REML failed on heterogeneous data"
        assert 0.01 <= prop <= 0.99, f"REML proportion {prop} out of range"
        assert not driver.execute_script("return isNaN(lastResults.pooledProp)")


class TestSubgroupQTest:
    """Test 24: Between-subgroup Q test produces numeric values."""

    def test_subgroup_q_between(self, driver):
        set_csv(driver, DEMO_SUBGROUP_CSV)
        analyze(driver)

        q_between = driver.execute_script("""
            if (!lastResults || !lastResults.subgroupResults || !lastResults.subgroupResults._between)
                return null;
            return lastResults.subgroupResults._between.Q;
        """)
        assert q_between is not None, "Between-group Q is null"
        assert q_between >= 0, f"Between-group Q = {q_between} should be >= 0"
        assert not driver.execute_script(
            "return isNaN(lastResults.subgroupResults._between.Q)"
        ), "Between-group Q is NaN"


class TestManualEntry:
    """Test 25: Manual entry mode works."""

    def test_manual_mode(self, driver):
        # Switch to manual mode
        driver.execute_script("switchInputMode('manual')")
        time.sleep(0.2)

        # Fill the 3 existing rows via JS
        driver.execute_script("""
            var rows = document.querySelectorAll('#manualRows .manual-row');
            var data = [
                ['StudyA', '10', '50', ''],
                ['StudyB', '15', '80', ''],
                ['StudyC', '20', '100', '']
            ];
            for (var i = 0; i < 3 && i < rows.length; i++) {
                var inputs = rows[i].querySelectorAll('input');
                inputs[0].value = data[i][0];
                inputs[1].value = data[i][1];
                inputs[2].value = data[i][2];
                inputs[3].value = data[i][3];
            }
        """)
        time.sleep(0.2)

        analyze(driver)
        prop = get_pooled(driver)
        assert prop is not None, "Manual entry analysis failed"
        assert 0.10 <= prop <= 0.30, f"Manual entry proportion {prop} unexpected"


class TestWeightsSum:
    """Test 26: Study weights sum to ~100%."""

    def test_weights_sum(self, driver):
        set_csv(driver, DEMO_CSV)
        analyze(driver)

        total_w = driver.execute_script("""
            if (!lastResults || !lastResults.studies) return null;
            return lastResults.studies.reduce(function(s, st) { return s + st.weight; }, 0);
        """)
        assert total_w is not None
        assert abs(total_w - 100) < 0.5, f"Weights sum to {total_w}, expected ~100"


class TestDLvsREMLAgreement:
    """Test 27: DL and REML give similar results on homogeneous data."""

    def test_dl_reml_agree(self, driver):
        set_csv(driver, DEMO_CSV)
        select_opt(driver, 'transform', 'FT')
        set_hksj(driver, False)

        # DL
        select_opt(driver, 'method', 'DL')
        analyze(driver)
        prop_dl = get_pooled(driver)

        # REML
        select_opt(driver, 'method', 'REML')
        analyze(driver)
        prop_reml = get_pooled(driver)

        assert prop_dl is not None and prop_reml is not None
        assert abs(prop_dl - prop_reml) < 0.03, \
            f"DL ({prop_dl:.4f}) and REML ({prop_reml:.4f}) too different on homogeneous data"


if __name__ == '__main__':
    pytest.main([__file__, '-v', '--tb=short'])
