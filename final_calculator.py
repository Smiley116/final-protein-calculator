import streamlit as st
import requests
import re
import time
import json
import base64
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# 工具函数
def clean_sequence(sequence):
    """清理序列"""
    cleaned = ''.join([char for char in sequence.upper() if char.isalpha()])
    valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
    
    # 只保留标准氨基酸
    cleaned = ''.join([aa for aa in cleaned if aa in valid_aa])
    
    return cleaned

def extract_sequence_from_input(input_str):
    """从输入文本中提取氨基酸序列或PDB ID"""
    if not input_str:
        return ""
    
    stripped_input = input_str.strip()
    
    # 检查是否为PDB ID (4个字符)
    if re.match(r'^[0-9a-zA-Z]{4}$', stripped_input):
        return stripped_input
    
    sequence_parts = []
    
    if input_str.startswith('>'):
        lines = input_str.split('\n')
        # 处理FASTA格式：跳过标题行，只处理后续行
        for line in lines[1:]:
            if line:
                sequence_parts.append(line)
    else:
        # 不是FASTA格式，直接处理整个字符串
        sequence_parts.append(input_str)
    
    sequence_content = ''.join(sequence_parts)
    return clean_sequence(sequence_content)

def analyze_sequence(sequence):
    prot = ProteinAnalysis(sequence)
    
    # 计算各项指标
    mw_da = prot.molecular_weight()
    mw_kda = mw_da / 1000
    pI = prot.isoelectric_point()
    
    # 消光系数计算
    # 获取原始消光系数 (M⁻¹cm⁻¹)
    ext_no_cys = prot.molar_extinction_coefficient()[0]  # 无二硫键
    ext_with_cys = prot.molar_extinction_coefficient()[1]  # 有二硫键
    
    # 计算Abs 0.1% (=1 g/L)的值
    # 正确公式：Abs 0.1% = 消光系数 / 分子量 × 1000
    abs_no_cys = ext_no_cys / mw_da  # 无二硫键
    abs_with_cys = ext_with_cys / mw_da # 有二硫键
    
    gravy = prot.gravy()
    aa_comp = prot.get_amino_acids_percent()
    
    return {
        'length': len(sequence),
        'mw_da': mw_da,
        'mw_kda': mw_kda,
        'pI': pI,
        'ext_no_cys': ext_no_cys,
        'ext_with_cys': ext_with_cys,
        'abs_no_cys': abs_no_cys,
        'abs_with_cys': abs_with_cys,
        'gravy': gravy,
        'aa_comp': aa_comp,
        'sequence': sequence
    }

def display_physicochemical_properties(result):
    """显示理化性质分析结果"""
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("### 基本信息")
        st.markdown(f"**序列长度**: {result['length']} 个氨基酸")
        st.markdown(f"**分子量**: {result['mw_kda']:.2f} kDa ({result['mw_da']:.0f} Da)")
        st.markdown(f"**等电点(pI)**: {result['pI']:.2f}")
        st.markdown(f"**平均疏水性(GRAVY)**: {result['gravy']:.3f}")
        
        # 解释疏水性
        if result['gravy'] > 0.5:
            hydrophobicity = "强疏水性"
        elif result['gravy'] > 0:
            hydrophobicity = "弱疏水性"
        elif result['gravy'] > -0.5:
            hydrophobicity = "弱亲水性"
        else:
            hydrophobicity = "强亲水性"
        st.markdown(f"**疏水性描述**: {hydrophobicity}")
    
    with col2:
        st.markdown("### 消光系数")
        # st.markdown(f"**无二硫键**: {result['ext_no_cys']:.0f} M⁻¹cm⁻¹")
        # st.markdown(f"**有二硫键**: {result['ext_with_cys']:.0f} M⁻¹cm⁻¹")
        st.markdown(f"**Abs 0.1% (1 mg/ml) - 无二硫键**: {result['abs_no_cys']:.3f}")
        st.markdown(f"**Abs 0.1% (1 mg/ml) - 有二硫键**: {result['abs_with_cys']:.3f}")
    
    # 氨基酸组成分析
    st.markdown("### 氨基酸组成分析")
    
    # 获取氨基酸组成并排序
    aa_comp = result['aa_comp']
    # 按照百分比从高到低排序
    sorted_aa = sorted(aa_comp.items(), key=lambda x: x[1], reverse=True)
    
    # 显示前10个最丰富的氨基酸
    st.markdown("#### 主要氨基酸（按丰度排序）")
    cols = st.columns(5)
    for j, (aa, percentage) in enumerate(sorted_aa[:10]):
        with cols[j % 5]:
            st.markdown(
                f"""
                <div class="metric-box">
                    <p style="font-size: 1.2rem; font-weight: bold; margin-bottom: 5px;">{aa}</p>
                    <p style="margin: 0; color: #333;">{percentage:.1%}</p>
                </div>
                """,
                unsafe_allow_html=True
            )
    
    # 显示疏水性和极性氨基酸统计
    hydrophobic = set(['A', 'I', 'L', 'M', 'F', 'W', 'V', 'P'])
    polar = set(['N', 'C', 'Q', 'S', 'T', 'Y'])
    charged = set(['R', 'H', 'K', 'D', 'E'])
    
    hydrophobic_count = sum(aa_comp.get(aa, 0) for aa in hydrophobic)
    polar_count = sum(aa_comp.get(aa, 0) for aa in polar)
    charged_count = sum(aa_comp.get(aa, 0) for aa in charged)
    
    st.markdown("#### 氨基酸分类统计")
    col1, col2, col3 = st.columns(3)
    with col1:
        st.markdown(
            f"""
            <div class="metric-box">
                <p style="margin: 0; color: #666;">疏水性氨基酸</p>
                <p style="font-size: 1.2rem; font-weight: bold; margin: 5px 0;">{hydrophobic_count:.1%}</p>
            </div>
            """,
            unsafe_allow_html=True
        )
    with col2:
        st.markdown(
            f"""
            <div class="metric-box">
                <p style="margin: 0; color: #666;">极性氨基酸</p>
                <p style="font-size: 1.2rem; font-weight: bold; margin: 5px 0;">{polar_count:.1%}</p>
            </div>
            """,
            unsafe_allow_html=True
        )
    with col3:
        st.markdown(
            f"""
            <div class="metric-box">
                <p style="margin: 0; color: #666;">带电氨基酸</p>
                <p style="font-size: 1.2rem; font-weight: bold; margin: 5px 0;">{charged_count:.1%}</p>
            </div>
            """,
            unsafe_allow_html=True
        )

# 设置页面配置
st.set_page_config(
    page_title="蛋白质理化性质分析工具",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# 自定义CSS样式
st.markdown("""
<style>
    .main-header {
        font-size: 2.5rem;
        color: #1f77b4;
        text-align: center;
        margin-bottom: 2rem;
    }
    .result-box {
        background-color: #f0f2f6;
        padding: 20px;
        border-radius: 10px;
        margin-top: 20px;
        border-left: 5px solid #1f77b4;
    }
    .metric-box {
        background-color: white;
        padding: 15px;
        border-radius: 10px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        margin: 5px;
    }
    .big-red-button {
        background-color: #e74c3c;
        color: white;
        font-size: 1.1rem;
        font-weight: bold;
        padding: 12px 24px;
        border-radius: 8px;
        border: none;
        width: 100%;
    }
    .big-red-button:hover {
        background-color: #c0392b;
        color: white;
    }
</style>
""", unsafe_allow_html=True)

# 标题
st.markdown('<div class="main-header">🧬 蛋白质理化性质分析工具</div>', unsafe_allow_html=True)

# Initialize session state variables for prediction
if 'prediction_status' not in st.session_state:
    st.session_state.prediction_status = []
if 'prediction_results' not in st.session_state:
    st.session_state.prediction_results = []
if 'api_settings' not in st.session_state:
    st.session_state.api_settings = {'use_api': False}
if 'last_analysis_result' not in st.session_state:
    st.session_state.last_analysis_result = None
if 'current_analysis_index' not in st.session_state:
    st.session_state.current_analysis_index = 0
if 'analysis_info' not in st.session_state:
    st.session_state.analysis_info = ""

# 侧边栏
with st.sidebar:
    st.header("📋 使用说明")
    st.markdown("""
    **输入方式：**
    - 🔍 **PDB ID**: 4位代码，如 `1crn`, `2abl`
    - 🧬 **氨基酸序列**: 直接粘贴序列
    - 📝 **FASTA格式**: 支持带>头的格式
    
    **分析内容：**
    - 分子量计算
    - 等电点(pI)
    - 消光系数
    - 疏水性(GRAVY)
    - 氨基酸组成分析
    - 🧬 蛋白质结构预测
    """)
    
    st.header("⚡ 示例")
    if st.button("加载胰岛素示例"):
        st.session_state.example_sequence = ">sp|P01308|INS_HUMAN Insulin\n" + \
            "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED\n" + \
            "LQVGQVELGGGPGAGSLQKRGIVEQCCTSICSLYQLENYCN"
    
    # API设置
    st.header("⚙️ API设置")
    st.checkbox("启用实际API", 
               value=st.session_state.api_settings['use_api'],
               disabled=True)
    st.info("API已预先配置，用户可直接使用蛋白质结构预测功能")
    st.markdown("""
    > **提示**: 管理员可修改代码中的API密钥配置
    """)
    
    st.header("🔗 关于")
    st.markdown("基于Python和BioPython开发的蛋白质分析与结构预测工具")

# 主界面 - 输入方式选择
input_method = st.radio(
    "选择输入方式:",
    ["直接输入氨基酸序列", "通过PDB ID分析"],
    horizontal=True,
    index=0  # 默认选择直接输入氨基酸序列
)

# 初始化序列列表
if 'sequences' not in st.session_state:
    st.session_state.sequences = [""]
# 确保预测状态和结果列表与序列列表长度一致 - 移至序列处理时统一处理

# 输入区域
if input_method == "通过PDB ID分析":
    col1, col2 = st.columns([1, 3])
    with col1:
        pdb_id = st.text_input("PDB ID:", placeholder="如: 1crn", key="pdb_input")
    with col2:
        if st.button("🔍 获取PDB示例", key="pdb_example"):
            pdb_id = "1crn"
            st.rerun()
else:
    st.subheader("🧬 氨基酸序列输入区")
    st.markdown("您可以添加多个氨基酸序列，系统将合并计算它们的总理化性质")
    
    # 初始化预测状态和结果
    if 'prediction_status' not in st.session_state:
        st.session_state.prediction_status = {}
    if 'prediction_results' not in st.session_state:
        st.session_state.prediction_results = {}
    # 初始化或更新API设置 - 强制使用真实API
    if 'api_settings' not in st.session_state:
        st.session_state.api_settings = {
            'api_key': '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',  # 使用用户提供的API key
            'use_api': True,  # 强制启用API，确保使用真实分析
            'api_url': 'https://health.api.nvidia.com/v1/biology/mit/boltz2/predict'  # NVIDIA的API URL
        }
    else:
        # 确保使用正确的设置，添加所有必需的键
        st.session_state.api_settings.setdefault('api_key', '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')
        st.session_state.api_settings['use_api'] = True  # 强制启用API
        st.session_state.api_settings['api_url'] = 'https://health.api.nvidia.com/v1/biology/mit/boltz2/predict'
        
    # 添加API状态显示（可选，用于调试）
    st.sidebar.markdown("### API设置状态")
    st.sidebar.markdown(f"- **API启用**: {st.session_state.api_settings['use_api']}")
    st.sidebar.markdown(f"- **API URL**: {st.session_state.api_settings['api_url']}")
    st.sidebar.markdown(f"- **API密钥**: {'已配置' if st.session_state.api_settings['api_key'] else '未配置'}")
    
    # 初始化分析相关状态
    if 'run_analysis' not in st.session_state:
        st.session_state.run_analysis = False
    if 'current_analysis_index' not in st.session_state:
        st.session_state.current_analysis_index = 0
    if 'last_analysis_result' not in st.session_state:
        st.session_state.last_analysis_result = None
    if 'analysis_info' not in st.session_state:
        st.session_state.analysis_info = ""
    
    # API设置已移除，默认使用API调用
    
    # 显示所有序列框
    for i in range(len(st.session_state.sequences)):
        # 序列标题
        st.markdown(f"### 🧬 序列 {i+1}")
        col1, col2 = st.columns([10, 1])
        with col1:
            st.session_state.sequences[i] = st.text_area(
                f"序列输入:",
                value=st.session_state.sequences[i],
                height=100,
                placeholder="请输入氨基酸序列或FASTA格式...",
                key=f"seq_input_{i}"
            )
        with col2:
            if i > 0:  # 不允许删除第一个序列框
                if st.button(f"🗑️", key=f"remove_seq_{i}"):
                    del st.session_state.sequences[i]
                    # 清理相关的预测状态
                    if i in st.session_state.prediction_status:
                        del st.session_state.prediction_status[i]
                    if i in st.session_state.prediction_results:
                        del st.session_state.prediction_results[i]
                    st.rerun()
    
        # 显示预测状态
        col_status = st.columns([1])[0]
        with col_status:
            if i in st.session_state.prediction_status:
                if st.session_state.prediction_status[i] == "running":
                    st.info("🔄 预测正在进行中...")
                elif st.session_state.prediction_status[i] == "completed":
                    st.success("✅ 预测完成!")
                elif st.session_state.prediction_status[i] == "error":
                    st.error("❌ 预测失败")
        
        # 大红色按钮行（只保留预测按钮样式）
        st.markdown(f"""<style>
            div[data-testid="stButton"]:has(button[data-testid="button-predict_seq_{i}"]) button {{
                background-color: #e74c3c;
                color: white;
                font-size: 1.1rem;
                font-weight: bold;
                padding: 12px 24px;
                border-radius: 8px;
                border: none;
                width: 100%;
            }}
            div[data-testid="stButton"]:has(button[data-testid="button-predict_seq_{i}"]) button:hover {{
                background-color: #c0392b;
                color: white;
            }}
        </style>""", unsafe_allow_html=True)
        
        # 结构预测按钮
        if st.button(f"🧬 预测结构", key=f"predict_seq_{i}"):
            original_sequence = st.session_state.sequences[i]
            sequence = extract_sequence_from_input(original_sequence)
            
            # 添加调试信息
            debug_info = []
            debug_info.append(f"原始序列 {i+1}: {original_sequence}")
            debug_info.append(f"清理后序列 {i+1}: {sequence}")
            debug_info.append(f"序列长度 {i+1}: {len(sequence)}")
            
            st.session_state.debug_info = debug_info
            
            if not sequence:
                st.warning(f"序列 {i+1}：请先输入有效的氨基酸序列")
            elif not re.match(r'^[0-9a-zA-Z]{4}$', sequence) and len(sequence) < 10:
                st.warning(f"序列 {i+1}：序列太短，至少需要10个氨基酸，当前长度为 {len(sequence)}")
            else:
                # Initialize prediction_status and prediction_results if not exists
                if 'prediction_status' not in st.session_state:
                    st.session_state.prediction_status = []
                if 'prediction_results' not in st.session_state:
                    st.session_state.prediction_results = []
                # 确保预测状态列表长度与序列列表一致
                while len(st.session_state.prediction_status) < len(st.session_state.sequences):
                    st.session_state.prediction_status.append("idle")
                # 确保预测结果列表长度与序列列表一致
                while len(st.session_state.prediction_results) < len(st.session_state.sequences):
                    st.session_state.prediction_results.append({})
                # 开始预测
                st.session_state.prediction_status[i] = "running"
                st.rerun()
        
        # 显示调试信息（如果有）
    if 'debug_info' in st.session_state and st.session_state.debug_info:
        with st.expander("调试信息", expanded=False):
            st.markdown("\n".join(st.session_state.debug_info))

    # 显示预测状态
    if 'prediction_status' in st.session_state:
        # 确保状态列表长度与序列列表一致
        while len(st.session_state.prediction_status) < len(st.session_state.sequences):
            st.session_state.prediction_status.append("idle")
        status = st.session_state.prediction_status[i]
        if status == "running":
            st.markdown("🔄 **预测进行中...**")
        elif status == "success":
            st.markdown("✅ **预测完成**")
        elif status == "error":
            st.markdown("❌ **预测失败**")
            
            if 'prediction_results' in st.session_state and i < len(st.session_state.prediction_results):
                result = st.session_state.prediction_results[i]
                if 'error' in result:
                    st.error(f"错误原因: {result['error']}")
        
        # 额外的错误信息检查，确保错误总是可见
        if status == "error" and 'prediction_results' in st.session_state and i < len(st.session_state.prediction_results):
            result = st.session_state.prediction_results[i]
            if 'error' in result:
                st.error(f"错误原因: {result['error']}")
    
    # 显示分析结果（如果有）
    if st.session_state.last_analysis_result is not None and st.session_state.get('current_analysis_index') == i:
        with st.expander(f"📊 {st.session_state.analysis_info}", expanded=True):
            result = st.session_state.last_analysis_result
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("### 基本信息")
                st.markdown(f"**序列长度**: {result['length']} 个氨基酸")
                st.markdown(f"**分子量**: {result['mw_kda']:.2f} kDa ({result['mw_da']:.0f} Da)")
                st.markdown(f"**等电点(pI)**: {result['pI']:.2f}")
                st.markdown(f"**平均疏水性(GRAVY)**: {result['gravy']:.3f}")
                
                # 解释疏水性
                if result['gravy'] > 0.5:
                    hydrophobicity = "强疏水性"
                elif result['gravy'] > 0:
                    hydrophobicity = "弱疏水性"
                elif result['gravy'] > -0.5:
                    hydrophobicity = "弱亲水性"
                else:
                    hydrophobicity = "强亲水性"
                st.markdown(f"**疏水性描述**: {hydrophobicity}")
            
            with col2:
                st.markdown("### 消光系数")
                st.markdown(f"**无二硫键**: {result['ext_no_cys']:.0f} M⁻¹cm⁻¹")
                st.markdown(f"**有二硫键**: {result['ext_with_cys']:.0f} M⁻¹cm⁻¹")
                st.markdown(f"**Abs 0.1% (1 mg/ml) - 无二硫键**: {result['abs_no_cys']:.3f}")
                st.markdown(f"**Abs 0.1% (1 mg/ml) - 有二硫键**: {result['abs_with_cys']:.3f}")
            
            # 氨基酸组成分析
            st.markdown("### 氨基酸组成分析")
            
            # 获取氨基酸组成并排序
            aa_comp = result['aa_comp']
            # 按照百分比从高到低排序
            sorted_aa = sorted(aa_comp.items(), key=lambda x: x[1], reverse=True)
            
            # 显示前10个最丰富的氨基酸
            st.markdown("#### 主要氨基酸（按丰度排序）")
            cols = st.columns(5)
            for j, (aa, percentage) in enumerate(sorted_aa[:10]):
                with cols[j % 5]:
                    st.markdown(
                        f"""
                        <div class="metric-box">
                            <p style="font-size: 1.2rem; font-weight: bold; margin-bottom: 5px;">{aa}</p>
                            <p style="margin: 0; color: #333;">{percentage:.1%}</p>
                        </div>
                        """,
                        unsafe_allow_html=True
                    )
            
            # 显示疏水性和极性氨基酸统计
            hydrophobic = set(['A', 'I', 'L', 'M', 'F', 'W', 'V', 'P'])
            polar = set(['N', 'C', 'Q', 'S', 'T', 'Y'])
            charged = set(['R', 'H', 'K', 'D', 'E'])
            
            hydrophobic_count = sum(aa_comp.get(aa, 0) for aa in hydrophobic)
            polar_count = sum(aa_comp.get(aa, 0) for aa in polar)
            charged_count = sum(aa_comp.get(aa, 0) for aa in charged)
            
            st.markdown("#### 氨基酸分类统计")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.markdown(
                    f"""
                    <div class="metric-box">
                        <p style="margin: 0; color: #666;">疏水性氨基酸</p>
                        <p style="font-size: 1.2rem; font-weight: bold; margin: 5px 0;">{hydrophobic_count:.1%}</p>
                    </div>
                    """,
                    unsafe_allow_html=True
                )
            with col2:
                st.markdown(
                    f"""
                    <div class="metric-box">
                        <p style="margin: 0; color: #666;">极性氨基酸</p>
                        <p style="font-size: 1.2rem; font-weight: bold; margin: 5px 0;">{polar_count:.1%}</p>
                    </div>
                    """,
                    unsafe_allow_html=True
                )
            with col3:
                st.markdown(
                    f"""
                    <div class="metric-box">
                        <p style="margin: 0; color: #666;">带电氨基酸</p>
                        <p style="font-size: 1.2rem; font-weight: bold; margin: 5px 0;">{charged_count:.1%}</p>
                    </div>
                    """,
                    unsafe_allow_html=True
                )
            
    st.markdown("---")
    
    # 显示预测结果
    if i in st.session_state.prediction_results:
        result = st.session_state.prediction_results[i]
        with st.expander("📊 预测结果详情", expanded=True):
                # 基本信息
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown(f"**预测置信度**: {result.get('confidence', 'N/A')}")
                    st.markdown(f"**预测时间**: {result.get('time', 'N/A')}")
                    if 'simulation' in result and result['simulation']:
                        st.info("📝 这是模拟数据，仅供演示使用")
                
                with col2:
                    # 显示结构质量指标（如果有）
                    if 'metrics' in result and result['metrics']:
                        st.markdown("**结构质量指标**:")
                        for metric_name, metric_value in result['metrics'].items():
                            metric_display = {
                                'plddt': 'pLDDT评分',
                                'tm_score': 'TM-Score',
                                'rmsd': 'RMSD (Å)'
                            }
                            st.markdown(f"- {metric_display.get(metric_name, metric_name)}: {metric_value}")
                
                # 提供结构文件下载
                if 'structure_data' in result and 'content' in result['structure_data'] and result['structure_data']['content']:
                    structure_content = result['structure_data']['content']
                    structure_format = result['structure_data'].get('format', 'pdb')
                    structure_id = result['structure_data'].get('structure_id', f'structure_{i+1}')
                    
                    # 根据格式设置文件扩展名和MIME类型
                    if structure_format.lower() == 'mmcif':
                        file_extension = 'cif'
                        mime_type = 'chemical/x-mmcif'
                        label = "💾 下载mmCIF文件"
                    else:  # 默认使用pdb
                        file_extension = 'pdb'
                        mime_type = 'chemical/x-pdb'
                        label = "💾 下载PDB文件"
                    
                    # 创建临时文件以下载
                    st.download_button(
                        label=label,
                        data=structure_content,
                        file_name=f"{structure_id}.{file_extension}",
                        mime=mime_type,
                        key=f"download_structure_{i}"
                    )
                    
                    # 简单的结构信息
                    atom_count = structure_content.count('ATOM')
                    st.markdown(f"**结构信息**: 包含 {atom_count} 个原子")
                    
                    # 提供结构可视化提示
                    st.markdown("""
                    **结构可视化提示**:
                    - 下载PDB文件后可使用PyMOL、UCSF、Chimera等工具查看
                    - 也可上传至 [RCSB 3D Viewer](https://rcsb.org/3d-view) 在线查看
                    """)
        
        st.markdown("---")
    
    # 添加序列按钮
    if st.button("➕ 添加序列", key="add_seq_main"):
        st.session_state.sequences.append("")
        st.rerun()
    
    # 操作按钮
    col1, col2, col3 = st.columns(3)
    with col1:
        if st.button("📝 加载示例序列", key="seq_example"):
            if not st.session_state.sequences[0]:  # 如果第一个序列框为空，则填充
                st.session_state.sequences[0] = ">sp|P01308|INS_HUMAN Insulin\n" + \
                    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED\n" + \
                    "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
            else:  # 否则添加到新的序列框
                st.session_state.sequences.append(">sp|P01308|INS_HUMAN Insulin\n" + \
                    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAED\n" + \
                    "LQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN")
            st.rerun()
    with col2:
        if st.button("🗑️ 清空所有序列", key="clear_all_seq"):
            st.session_state.sequences = [""]
            st.rerun()
    with col3:
        st.markdown("**相关资源：** [NVIDIA Boltz-2 API](https://build.nvidia.com/mit/boltz2/apireference)")
        if st.button("🔄 预测所有序列", key="predict_all_seq"):
            # 验证所有序列
            all_valid = True
            debug_info = []
            for i, seq_input in enumerate(st.session_state.sequences):
                sequence = extract_sequence_from_input(seq_input)
                debug_info.append(f"原始输入序列 {i+1}: {seq_input}")
                debug_info.append(f"提取后序列 {i+1}: {sequence}")
                debug_info.append(f"提取后长度 {i+1}: {len(sequence)}")
                if not sequence:
                    st.warning(f"序列 {i+1} 为空，请输入有效序列")
                    all_valid = False
                elif not re.match(r'^[0-9a-zA-Z]{4}$', sequence) and len(sequence) < 10:
                    st.warning(f"序列 {i+1} 太短，至少需要10个氨基酸，当前长度为 {len(sequence)}")
                    all_valid = False
            
            # 保存调试信息
            st.session_state.debug_info = debug_info
            
            if all_valid:
                # 确保预测状态列表长度与序列列表一致
                while len(st.session_state.prediction_status) < len(st.session_state.sequences):
                    st.session_state.prediction_status.append("idle")
                # 设置所有序列的预测状态为运行中
                for i in range(len(st.session_state.sequences)):
                    st.session_state.prediction_status[i] = "running"
                st.rerun()

# 示例序列处理（保持兼容性）
if 'example_sequence' in st.session_state:
    st.session_state.sequences[0] = st.session_state.example_sequence
    del st.session_state.example_sequence  # 删除临时存储

# 添加上方统一分析按钮
st.markdown("""
    <style>
        div[data-testid="stButton"]:has(button[data-testid="button-analyze_all"]) button {{
            background-color: #1f77b4;
            color: white;
            font-size: 1.1rem;
            font-weight: bold;
            padding: 12px 24px;
            border-radius: 8px;
            border: none;
            width: 100%;
            margin-top: 20px;
        }}
        div[data-testid="stButton"]:has(button[data-testid="button-analyze_all"]) button:hover {{
            background-color: #1565c0;
            color: white;
        }}
    </style>
""", unsafe_allow_html=True)

if st.button("🧪 分析所有输入序列的理化性质", key="analyze_all_bottom"):
    merged_sequence = ""
    valid_sequences_found = False
    individual_results = []  # 存储单独结果

    for i in range(len(st.session_state.sequences)):
        sequence = extract_sequence_from_input(st.session_state.sequences[i])
        if sequence:
            valid_sequences_found = True
            merged_sequence += sequence  # 合并所有有效序列
            
            # 单独计算并保存结果
            single_result = analyze_sequence(sequence)
            individual_results.append({
                'index': i,
                'result': single_result
            })
            
    if not valid_sequences_found:
        st.warning("请先输入有效的氨基酸序列")
    else:
        # 使用合并后的序列进行分析
        st.session_state.current_analysis_index = -1  # 设置为-1表示全局合并分析
        st.session_state.analysis_info = "合并所有序列的理化性质分析结果"
        # 调用分析函数
        result = analyze_sequence(merged_sequence)
        st.session_state.last_analysis_result = result
        # 保存单独分析结果
        st.session_state.individual_analysis_results = individual_results

# 显示合并分析结果
if st.session_state.last_analysis_result is not None and st.session_state.get('current_analysis_index') == -1:
    # 显示合并结果
    with st.expander(f"📊 {st.session_state.analysis_info}", expanded=True):
        display_physicochemical_properties(st.session_state.last_analysis_result)
        
    # 显示单独结果（如果存在）
    if 'individual_analysis_results' in st.session_state and st.session_state.individual_analysis_results:
        st.markdown("### 🧬 各序列单独理化性质")
        for item in st.session_state.individual_analysis_results:
            idx = item['index']
            res = item['result']
            # 使用 expander 默认折叠，实现用户要求的“点击可以显示”
            with st.expander(f"序列 {idx+1} 详情", expanded=False):
                display_physicochemical_properties(res)

# 不再使用全局分析按钮，已移至每个序列的独立按钮

# 删除全局分析按钮的逻辑
# if analyze_btn:
#    相关代码已移至每个序列的独立按钮中

def is_valid_pdb_id(pdb_id):
    """验证PDB ID格式"""
    return len(pdb_id) == 4 and re.match(r'^[0-9a-zA-Z]{4}$', pdb_id)

def get_sequence_from_pdb(pdb_id):
    """从PDB获取序列"""
    if not is_valid_pdb_id(pdb_id):
        st.error(f"{pdb_id} 不是有效的PDB ID")
        return None
    
    fasta_url = f"https://www.rcsb.org/fasta/entry/{pdb_id.upper()}/download"
    
    try:
        response = requests.get(fasta_url, timeout=10)
        response.raise_for_status()
        
        lines = [line.strip() for line in response.text.split('\n') if line.strip()]
        sequence_lines = [line for line in lines[1:] if "|" not in line]  
        raw_sequence = ''.join(sequence_lines).replace(' ', '').replace('\n', '')
        
        return clean_sequence(raw_sequence)
        
    except requests.exceptions.RequestException as e:
        st.error(f"PDB序列下载失败: {e}")
        return None

def mock_protein_structure_prediction(sequence):
    """
    模拟蛋白质结构预测函数
    用于在没有实际API密钥时演示功能
    
    Args:
        sequence: 氨基酸序列
    
    Returns:
        dict: 包含模拟预测结果的字典
    """
    # 模拟预测延迟
    time.sleep(2)
    
    # 模拟预测结果
    import random
    confidence = round(random.uniform(70, 99), 1)
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    
    # 生成模拟的PDB格式结构数据
    mock_pdb_data = """HEADER    SIMULATED PROTEIN STRUCTURE    {timestamp}
TITLE     MOCK PREDICTION RESULT
COMPND    MOCK PROTEIN
SOURCE    SIMULATED BY TRAE AI
KEYWDS    MOCK, STRUCTURE, PREDICTION
EXPDTA    MOCK DATA
AUTHOR    TRAE AI
REMARK    1 AUTH GENERATED BY MOCK PROTEIN STRUCTURE PREDICTION
ATOM      1  N   GLY A   1      12.431  -1.025   0.762  1.00 99.99           N
ATOM      2  CA  GLY A   1      11.766   0.098   0.228  1.00 99.99           C
ATOM      3  C   GLY A   1      12.374   1.319   0.663  1.00 99.99           C
ATOM      4  O   GLY A   1      13.552   1.417   0.385  1.00 99.99           O
ATOM      5  N   ALA A   2      11.562   2.285   1.334  1.00 99.99           N
ATOM      6  CA  ALA A   2      12.043   3.560   1.869  1.00 99.99           C
ATOM      7  C   ALA A   2      11.111   4.638   1.661  1.00 99.99           C
ATOM      8  O   ALA A   2      10.852   5.738   2.238  1.00 99.99           O
ATOM      9  CB  ALA A   2      13.563   3.810   1.526  1.00 99.99           C
ATOM     10  N   SER A   3      10.541   4.408   0.474  1.00 99.99           N
ATOM     11  CA  SER A   3       9.658   5.368   0.145  1.00 99.99           C
ATOM     12  C   SER A   3       8.264   4.957  -0.115  1.00 99.99           C
ATOM     13  O   SER A   3       7.228   5.592  -0.293  1.00 99.99           O
ATOM     14  CB  SER A   3      10.023   6.731  -0.365  1.00 99.99           C
ATOM     15  OG  SER A   3      11.284   7.120  -0.705  1.00 99.99           O
ENDMDL
""".format(timestamp=timestamp)
    
    # 生成模拟的结构质量指标
    mock_metrics = {
        "plddt": round(random.uniform(70, 95), 1),
        "tm_score": round(random.uniform(0.7, 0.95), 3),
        "rmsd": round(random.uniform(0.5, 3.0), 2)
    }
    
    return {
        "confidence": f"{confidence}%",
        "time": timestamp,
        "metrics": mock_metrics,
        "structure_data": {
            "content": mock_pdb_data,
            "structure_id": f"mock_{time.strftime('%Y%m%d_%H%M%S')}"
        },
        "message": "这是一个模拟的蛋白质结构预测结果，用于演示功能",
        "simulation": True
    }

def mock_calculate_affinity(sequence1, sequence2):
    """
    模拟蛋白质亲和度计算函数
    用于在没有实际Boltz API密钥时演示功能
    
    Args:
        sequence1: 第一个蛋白质氨基酸序列
        sequence2: 第二个蛋白质氨基酸序列
    
    Returns:
        dict: 包含模拟亲和度结果的字典
    """
    # 模拟API延迟
    time.sleep(1.5)
    
    import random
    
    # 生成模拟亲和度数据
    affinity_score = round(random.uniform(0.5, 1.0), 4)
    binding_energy = round(random.uniform(-15, -1), 2)
    dissociation_constant = round(random.uniform(1e-12, 1e-6), 12)
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    
    # 生成结合强度描述
    if affinity_score > 0.85:
        binding_strength = "极强"
    elif affinity_score > 0.7:
        binding_strength = "强"
    elif affinity_score > 0.55:
        binding_strength = "中等"
    else:
        binding_strength = "弱"
    
    return {
        "success": True,
        "sequence1": sequence1,
        "sequence2": sequence2,
        "affinity_score": affinity_score,
        "binding_energy": binding_energy,  # 结合能 (kcal/mol)
        "kd_value": dissociation_constant,  # 解离常数 (M)
        "binding_strength": binding_strength,
        "time": timestamp,
        "simulation": True,
        "method": "mock_binding_prediction"
    }

def api_protein_structure_prediction(sequence, api_key, api_url):
    """
    使用API进行蛋白质结构预测
    
    Args:
        sequence: 氨基酸序列
        api_key: API密钥
        api_url: API端点URL
    
    Returns:
        dict: 包含预测结果的字典
    
    Raises:
        Exception: 当API请求失败时
    """
    # 初始化调试信息列表
    debug_info = []
    debug_info.append("🔍 开始预测调试信息:")
    debug_info.append(f"- 序列长度: {len(sequence)}")
    debug_info.append(f"- API URL: {api_url}")
    debug_info.append(f"- API Key格式检查: {'有效' if api_key.startswith('nvapi-') else '无效'}")
    
    # 构建请求头 - 根据示例代码添加必要的轮询参数
    NVCF_POLL_SECONDS = 300
    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}",
        "NVCF-POLL-SECONDS": str(NVCF_POLL_SECONDS)
    }
    
    # 构建请求体 - 根据示例代码的正确格式
    payload = {
        "polymers": [
            {
                "id": "A",
                "molecule_type": "protein",
                "sequence": sequence
            }
        ],
        "recycling_steps": 1,
        "sampling_steps": 50,
        "diffusion_samples": 3,
        "step_scale": 1.2,
        "without_potentials": True
    }
    try:
        # 发送预测请求 - 增加超时时间并添加更详细的调试信息
        debug_info.append(f"- 正在连接到API: {api_url}")
        debug_info.append("- 正在发送请求数据，请稍候...")
        debug_info.append("- 注意：NVIDIA API可能需要较长响应时间(1-5分钟)")
        
        # 增加超时时间到300秒以匹配NVCF-POLL-SECONDS设置，应对API可能的长时间响应
        response = requests.post(api_url, json=payload, headers=headers, timeout=300)
        
        # 添加更多调试信息
        debug_info.append(f"- HTTP状态码: {response.status_code}")
        debug_info.append(f"- 响应头: {dict(response.headers)}")
        
        # 保存调试信息到session_state
        if 'debug_info' not in st.session_state:
            st.session_state.debug_info = []
        st.session_state.debug_info = debug_info
        
        # 处理202 Accepted响应 - 根据示例代码实现轮询逻辑
        if response.status_code == 202:
            debug_info.append("- 收到202 Accepted响应，开始轮询任务状态")
            task_id = response.headers.get("nvcf-reqid")
            
            if not task_id:
                raise Exception("未从202响应中获取到task_id")
                
            debug_info.append(f"- 获取到task_id: {task_id}")
            status_url = f"https://api.nvcf.nvidia.com/v2/nvcf/pexec/status/{task_id}"
            debug_info.append(f"- 状态查询URL: {status_url}")
            
            # 轮询状态
            max_retries = 30
            retry_count = 0
            while retry_count < max_retries:
                retry_count += 1
                debug_info.append(f"- 轮询尝试 {retry_count}/{max_retries}")
                
                status_response = requests.get(status_url, headers=headers, timeout=120)
                debug_info.append(f"  - 状态响应码: {status_response.status_code}")
                
                if status_response.status_code == 200:
                    debug_info.append("  - 任务完成，获取到结果")
                    # 更新调试信息
                    st.session_state.debug_info = debug_info
                    # 使用状态响应继续处理
                    response = status_response
                    break
                elif status_response.status_code in [400, 401, 404, 422, 500]:
                    error_msg = f"轮询任务状态失败 (状态码: {status_response.status_code})\n{status_response.text}"
                    debug_info.append(f"  - 轮询失败: {error_msg}")
                    st.session_state.debug_info = debug_info
                    raise Exception(error_msg)
                
                # 等待后重试
                wait_time = 20
                debug_info.append(f"  - 任务尚未完成，{wait_time}秒后重试")
                time.sleep(wait_time)
            
            if retry_count >= max_retries:
                raise Exception(f"轮询超时，已尝试 {max_retries} 次")
        
        # 检查HTTP响应状态
        if response.status_code != 200:
            error_msg = f"API请求失败 (状态码: {response.status_code})\n"
            try:
                error_data = response.json()
                debug_info.append(f"- 错误响应JSON: {error_data}")
                if "error" in error_data:
                    error_msg += f"错误详情: {error_data['error']}"
                elif "message" in error_data:
                    error_msg += f"错误信息: {error_data['message']}"
                else:
                    error_msg += f"响应内容: {response.text[:500]}..."
            except Exception as json_error:
                debug_info.append(f"- JSON解析错误: {str(json_error)}")
                error_msg += f"响应内容: {response.text[:500]}..."
            
            # 更新调试信息
            st.session_state.debug_info = debug_info
            raise Exception(error_msg)
        
        # 解析响应
        result = response.json()
        debug_info.append("- 响应JSON格式正确")
        debug_info.append(f"- 响应内容概览: {list(result.keys())}")
        
        # 添加更详细的响应结构调试信息
        debug_info.append("- 响应详细结构:")
        for key, value in result.items():
            if isinstance(value, dict):
                debug_info.append(f"  * {key}: {list(value.keys())}")
            elif isinstance(value, list):
                debug_info.append(f"  * {key}: 列表，长度={len(value)}")
                # 如果是structures列表，显示第一个元素的信息
                if key == 'structures' and value:
                    first_item = value[0]
                    if isinstance(first_item, dict):
                        debug_info.append(f"    - 第一个结构: {list(first_item.keys())}")
            else:
                debug_info.append(f"  * {key}: {type(value).__name__}")
        
        # 更新调试信息
        st.session_state.debug_info = debug_info
        
        # 标准化结果格式 - 根据示例代码中的响应结构
        structure_content = ""
        structure_format = "mmcif"  # 默认格式
        confidence_value = "未知"
        
        # 1. 从structures列表中提取结构数据（根据示例代码）
        if 'structures' in result and isinstance(result['structures'], list) and len(result['structures']) > 0:
            debug_info.append(f"- 发现structures列表，包含 {len(result['structures'])} 个结构")
            # 获取第一个结构
            first_structure = result['structures'][0]
            
            if isinstance(first_structure, dict):
                # 根据示例代码，结构数据存储在'structure'字段中
                if 'structure' in first_structure:
                    structure_content = first_structure['structure']
                    debug_info.append("- 从'structure'字段提取结构数据")
                    
                    # 确定格式
                    if 'format' in first_structure:
                        structure_format = first_structure['format']
                        debug_info.append(f"- 格式从'format'字段确定: {structure_format}")
                    else:
                        # 尝试根据内容判断
                        if structure_content.strip().startswith('HEADER'):
                            structure_format = "pdb"
                        else:
                            structure_format = "mmcif"
                        debug_info.append(f"- 自动判断格式: {structure_format}")
                else:
                    debug_info.append(f"- 结构字典中没有'structure'字段，键列表: {list(first_structure.keys())}")
                    # 尝试其他可能的字段
                    for key in ['content', 'pdb', 'mmcif']:
                        if key in first_structure:
                            structure_content = first_structure[key]
                            structure_format = key if key in ['pdb', 'mmcif'] else 'unknown'
                            debug_info.append(f"- 从'{key}'字段提取结构数据")
                            break
            else:
                debug_info.append(f"- 结构项不是字典，类型: {type(first_structure).__name__}")
        else:
            debug_info.append("- 未发现有效的structures列表")
            # 尝试其他可能的位置
            for path in ['prediction.structure', 'prediction', '']:
                current = result
                if path:
                    parts = path.split('.')
                    for part in parts:
                        if isinstance(current, dict) and part in current:
                            current = current[part]
                        else:
                            current = None
                            break
                
                if current:
                    if isinstance(current, dict):
                        for key in ['structure', 'content', 'pdb', 'mmcif']:
                            if key in current:
                                structure_content = current[key]
                                structure_format = key if key in ['pdb', 'mmcif'] else 'unknown'
                                debug_info.append(f"- 从{path}.{key}提取结构数据")
                                break
                    elif isinstance(current, str):
                        structure_content = current
                        debug_info.append(f"- 从{path}提取结构字符串")
        
        # 2. 提取置信度信息（根据示例代码）
        if 'confidence_scores' in result and isinstance(result['confidence_scores'], list) and len(result['confidence_scores']) > 0:
            confidence_value = f"{result['confidence_scores'][0]:.2f}"
            debug_info.append(f"- 从confidence_scores提取置信度: {confidence_value}")
        else:
            debug_info.append("- 未找到confidence_scores字段")
            # 尝试其他可能的置信度字段
            for score_key in ['iptm_scores', 'ptm_scores', 'confidence']:
                if score_key in result:
                    if isinstance(result[score_key], list) and result[score_key]:
                        confidence_value = f"{result[score_key][0]:.2f}"
                        debug_info.append(f"- 从{score_key}提取置信度: {confidence_value}")
                    elif isinstance(result[score_key], (int, float)):
                        confidence_value = f"{result[score_key]:.2f}"
                        debug_info.append(f"- 从{score_key}提取置信度: {confidence_value}")
                    break
        
        # 添加格式信息到调试日志
        debug_info.append(f"- 最终检测到的结构格式: {structure_format}")
        debug_info.append(f"- 结构数据长度: {len(structure_content) if structure_content else 0} 字符")
        
        # 提取metrics信息
        metrics = {}
        # 尝试从不同位置提取metrics
        if 'metrics' in result:
            metrics = result['metrics']
            debug_info.append("- 从根级提取metrics")
        elif 'prediction' in result and isinstance(result['prediction'], dict) and 'metrics' in result['prediction']:
            metrics = result['prediction']['metrics']
            debug_info.append("- 从prediction提取metrics")
        
        # 添加metrics信息到调试日志
        if metrics:
            debug_info.append(f"- 提取到metrics: {list(metrics.keys())}")
        else:
            debug_info.append("- 未提取到metrics")
        
        # 标准化结果，使用新提取的置信度值和结构数据
        # 保存更多原始响应数据用于调试
        raw_response_data = {
            "has_structures": 'structures' in result and isinstance(result['structures'], list),
            "response_keys": list(result.keys()),
            "has_confidence_scores": 'confidence_scores' in result and isinstance(result['confidence_scores'], list),
            "has_metrics": 'metrics' in result and isinstance(result['metrics'], dict)
        }
        
        # 添加structures列表的前几个元素（如果存在）
        if 'structures' in result and isinstance(result['structures'], list):
            raw_response_data['structures_count'] = len(result['structures'])
            # 保存第一个结构的关键信息
            if result['structures'] and isinstance(result['structures'][0], dict):
                first_struct = result['structures'][0]
                raw_response_data['first_structure_keys'] = list(first_struct.keys())
                # 保存格式信息
                if 'format' in first_struct:
                    raw_response_data['first_structure_format'] = first_struct['format']
        
        # 添加置信度信息
        if 'confidence_scores' in result:
            raw_response_data['confidence_scores_type'] = type(result['confidence_scores']).__name__
            if isinstance(result['confidence_scores'], list):
                raw_response_data['confidence_scores_count'] = len(result['confidence_scores'])
        
        # 添加metrics信息
        if 'metrics' in result:
            raw_response_data['metrics_keys'] = list(result['metrics'].keys()) if isinstance(result['metrics'], dict) else None
        
        standardized_result = {
            "confidence": confidence_value,
            "time": time.strftime("%Y-%m-%d %H:%M:%S"),
            "structure_content": structure_content,
            "structure_format": structure_format,
            "metrics": metrics,
            "structure_data": {
                "content": structure_content,
                "format": structure_format,
                "structure_id": f"nvidia_{time.strftime('%Y%m%d_%H%M%S')}"
            },
            "simulation": False,
            "raw_response": raw_response_data  # 保存更详细的原始响应信息
        }
        
        # 添加标准化结果的调试信息
        debug_info.append("- 标准化结果概览:")
        debug_info.append(f"  * 置信度: {standardized_result['confidence']}")
        debug_info.append(f"  * 指标数量: {len(standardized_result['metrics'])}")
        debug_info.append(f"  * 结构格式: {standardized_result['structure_data']['format']}")
        debug_info.append(f"  * 结构ID: {standardized_result['structure_data']['structure_id']}")
        st.session_state.debug_info = debug_info
        
        return standardized_result
        
    except requests.exceptions.Timeout:
        debug_info.append("- 错误类型: 请求超时错误")
        st.session_state.debug_info = debug_info
        raise Exception("API请求超时，请稍后再试")
    except requests.exceptions.ConnectionError:
        debug_info.append("- 错误类型: 连接错误")
        st.session_state.debug_info = debug_info
        raise Exception("无法连接到API服务器，请检查网络连接")
    except requests.exceptions.RequestException as e:
        debug_info.append(f"- 错误类型: 请求异常")
        debug_info.append(f"- 错误详情: {str(e)}")
        st.session_state.debug_info = debug_info
        raise Exception(f"API请求错误: {str(e)}")
    except Exception as e:
        debug_info.append(f"- 错误类型: 其他异常")
        debug_info.append(f"- 错误详情: {str(e)}")
        st.session_state.debug_info = debug_info
        raise Exception(f"预测过程中出错: {str(e)}")

# 处理正在进行的预测任务
def process_pending_predictions():
    """处理所有正在进行的预测任务"""
    # 确保序列列表存在
    if 'sequences' not in st.session_state:
        st.warning("序列列表不存在")
        return False
    
    if 'prediction_status' in st.session_state:
        # 初始化预测结果列表
        if 'prediction_results' not in st.session_state:
            st.session_state.prediction_results = []
        
        # 确保prediction_status和sequences长度一致
        while len(st.session_state.prediction_status) < len(st.session_state.sequences):
            st.session_state.prediction_status.append("idle")
        
        # 标记是否需要重新运行
        need_rerun = False
        
        # 添加调试信息
        debug_info = []
        debug_info.append(f"开始处理预测任务 - 总共 {len(st.session_state.prediction_status)} 个任务")
        debug_info.append(f"序列列表长度: {len(st.session_state.sequences)}")
        
        # 初始化进度列表
        progress = [0] * len(st.session_state.prediction_status)
        
        # 遍历所有状态为running的任务
        for seq_idx, status in enumerate(st.session_state.prediction_status):
            debug_info.append(f"任务 {seq_idx+1} 状态: {status}")
            if status == "running" and seq_idx < len(st.session_state.sequences):
                try:
                    # 获取序列并进行预测
                    original_sequence = st.session_state.sequences[seq_idx]
                    debug_info.append(f"任务 {seq_idx+1} - 原始序列: {original_sequence}")
                    
                    sequence = extract_sequence_from_input(original_sequence)
                    debug_info.append(f"任务 {seq_idx+1} - 提取后序列: {sequence}")
                    debug_info.append(f"任务 {seq_idx+1} - 序列长度: {len(sequence) if sequence else 0}")
                    
                    # 首先检查序列是否有效
                    if not sequence:
                        # 序列无效，取消预测
                        st.session_state.prediction_status[seq_idx] = "error"
                        error_msg = f"序列 {seq_idx+1} 无效或为空"
                        st.error(error_msg)
                        debug_info.append(f"任务 {seq_idx+1} - 序列无效，标记为错误")
                        
                        # 确保预测结果列表长度足够
                        while len(st.session_state.prediction_results) <= seq_idx:
                            st.session_state.prediction_results.append({})
                        
                        # 保存错误信息到结果中
                        st.session_state.prediction_results[seq_idx] = {
                            'error': '无效或空序列',
                            'error_type': 'invalid_sequence',
                            'sequence_index': seq_idx
                        }
                        
                        # 更新进度
                        progress[seq_idx] = 100
                        need_rerun = True  # 需要重新运行以更新UI
                        continue  # 继续处理下一个任务
                    
                    # 序列有效，继续处理
                    sequence_to_predict = sequence
                    
                    # 检查是否是PDB ID
                    if re.match(r'^[0-9a-zA-Z]{4}$', sequence):
                        # 从PDB数据库获取实际序列
                        try:
                            pdb_sequence = get_sequence_from_pdb(sequence)
                            if pdb_sequence:
                                sequence_to_predict = pdb_sequence
                                debug_info.append(f"任务 {seq_idx+1} - 从PDB获取序列成功，长度: {len(pdb_sequence)}")
                            else:
                                # 获取序列失败
                                st.session_state.prediction_status[seq_idx] = "error"
                                error_msg = f"无法获取PDB ID {sequence} 的序列"
                                st.error(error_msg)
                                debug_info.append(f"任务 {seq_idx+1} - 获取PDB序列失败")
                                
                                # 确保预测结果列表长度足够
                                while len(st.session_state.prediction_results) <= seq_idx:
                                    st.session_state.prediction_results.append({})
                                
                                # 保存错误信息到结果中
                                st.session_state.prediction_results[seq_idx] = {
                                    'error': f'无法获取PDB ID {sequence} 的序列',
                                    'error_type': 'pdb_sequence_error',
                                    'sequence_index': seq_idx
                                }
                                
                                need_rerun = True
                                continue  # 继续处理下一个任务
                        except Exception as pdb_error:
                            st.session_state.prediction_status[seq_idx] = "error"
                            error_msg = f"获取PDB序列时出错: {str(pdb_error)}"
                            st.error(error_msg)
                            debug_info.append(f"任务 {seq_idx+1} - PDB错误: {str(pdb_error)}")
                            
                            # 确保预测结果列表长度足够
                            while len(st.session_state.prediction_results) <= seq_idx:
                                st.session_state.prediction_results.append({})
                            
                            # 保存错误信息到结果中
                            st.session_state.prediction_results[seq_idx] = {
                                'error': f'获取PDB序列时出错: {str(pdb_error)}',
                                'error_type': 'pdb_exception',
                                'sequence_index': seq_idx
                            }
                            
                            need_rerun = True
                            continue
                    
                    # 不是PDB ID，检查序列长度
                    if len(sequence_to_predict) < 10:
                        st.session_state.prediction_status[seq_idx] = "error"
                        error_msg = f"序列 {seq_idx+1} 太短，至少需要10个氨基酸"
                        st.error(error_msg)
                        debug_info.append(f"任务 {seq_idx+1} - 序列太短，跳过")
                        
                        # 确保预测结果列表长度足够
                        while len(st.session_state.prediction_results) <= seq_idx:
                            st.session_state.prediction_results.append({})
                        
                        # 保存错误信息到结果中
                        st.session_state.prediction_results[seq_idx] = {
                            'error': f'序列太短，至少需要10个氨基酸，当前长度为 {len(sequence_to_predict)}',
                            'error_type': 'sequence_too_short',
                            'sequence': sequence_to_predict,
                            'sequence_index': seq_idx
                        }
                        
                        need_rerun = True
                        continue  # 继续处理下一个任务
                    
                    # 确保api_settings已初始化
                    if 'api_settings' not in st.session_state:
                        st.session_state.api_settings = {
                            'use_api': False,
                            'api_key': '',
                            'api_url': 'https://health.api.nvidia.com/v1/biology/mit/boltz2/predict'
                        }
                    
                    debug_info.append(f"任务 {seq_idx+1} - API使用状态: {st.session_state.api_settings['use_api']}")
                    
                    # 执行预测
                    try:
                        # 预测处理
                        if st.session_state.api_settings['use_api'] and st.session_state.api_settings['api_key']:
                            # 使用实际API进行预测
                            debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 调用实际API进行预测")
                            result = api_protein_structure_prediction(
                                sequence_to_predict,
                                st.session_state.api_settings['api_key'],
                                st.session_state.api_settings['api_url']
                            )
                        else:
                            # 使用模拟函数进行预测
                            debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 使用模拟数据进行预测")
                            # 添加详细的调用信息
                            debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 调用mock_protein_structure_prediction")
                            debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 输入序列长度: {len(sequence_to_predict)}")
                            result = mock_protein_structure_prediction(sequence_to_predict)
                        
                        # 详细检查预测结果
                        debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 预测完成，结果类型: {type(result)}")
                        
                        if not isinstance(result, dict):
                            raise TypeError(f"预测结果应该是字典类型，但实际是 {type(result).__name__}")
                        
                        # 检查结果中的必要键
                        debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 结果包含的键: {list(result.keys())}")
                        
                        # 检查是否有structure_data键
                        if 'structure_data' not in result:
                            debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 警告: 结果中缺少'structure_data'键")
                            # 尝试查找可能的替代键
                            structure_keys = [k for k in result.keys() if 'structure' in k.lower() or 'pdb' in k.lower()]
                            debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 可能的结构数据键: {structure_keys}")
                        else:
                            # 检查structure_data的结构
                            structure_data = result['structure_data']
                            debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - structure_data类型: {type(structure_data)}")
                            
                            if isinstance(structure_data, dict):
                                debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - structure_data中的键: {list(structure_data.keys())}")
                                # 检查是否有content或pdb键
                                if 'content' not in structure_data and 'pdb' not in structure_data:
                                    debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 警告: structure_data中缺少'content'或'pdb'键")
                        
                        # 添加序列信息到结果
                        result['sequence'] = sequence_to_predict
                        result['sequence_index'] = seq_idx
                        result['original_sequence'] = original_sequence
                        result['prediction_timestamp'] = time.time()
                        
                        # 确保预测结果列表长度足够
                        while len(st.session_state.prediction_results) <= seq_idx:
                            st.session_state.prediction_results.append({})
                        
                        # 保存预测结果
                        st.session_state.prediction_results[seq_idx] = result
                        st.session_state.prediction_status[seq_idx] = "success"
                        debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 预测成功，结果已保存")
                        need_rerun = True
                        
                        # 更新进度
                        progress[seq_idx] = 100
                    except Exception as prediction_error:
                        # 处理预测错误，收集详细的异常信息
                        import traceback
                        error_type = type(prediction_error).__name__
                        error_message = str(prediction_error)
                        error_trace = traceback.format_exc()
                        
                        st.session_state.prediction_status[seq_idx] = "error"
                        error_msg = f"预测失败: {error_message}"
                        st.error(error_msg)
                        
                        # 记录详细的错误信息
                        debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 预测错误")
                        debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 错误类型: {error_type}")
                        debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 错误消息: {error_message}")
                        debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 错误堆栈: {error_trace[:500]}..." if len(error_trace) > 500 else f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 错误堆栈: {error_trace}")
                        
                        # 确保预测结果列表长度足够
                        while len(st.session_state.prediction_results) <= seq_idx:
                            st.session_state.prediction_results.append({})
                        
                        # 保存详细的错误信息到结果中
                        st.session_state.prediction_results[seq_idx] = {
                            'error': f'预测失败: {error_message}',
                            'error_type': error_type,
                            'error_trace': error_trace,
                            'sequence_index': seq_idx,
                            'sequence': sequence_to_predict,
                            'error_timestamp': time.time()
                        }
                        
                        need_rerun = True
                        debug_info.append(f"[{time.strftime('%H:%M:%S')}] 任务 {seq_idx+1} - 错误信息已保存到索引 {seq_idx}")
                        
                        # 更新进度
                        progress[seq_idx] = 100
                    
                except Exception as task_error:
                    # 捕获任务级别的异常
                    st.session_state.prediction_status[seq_idx] = "error"
                    error_msg = f"处理任务时出错: {str(task_error)}"
                    st.error(error_msg)
                    debug_info.append(f"任务 {seq_idx+1} - 任务级错误: {str(task_error)}")
                    
                    # 确保预测结果列表长度足够
                    while len(st.session_state.prediction_results) <= seq_idx:
                        st.session_state.prediction_results.append({})
                    
                    # 保存错误信息到结果中
                    st.session_state.prediction_results[seq_idx] = {
                        'error': f'任务处理错误: {str(task_error)}',
                        'error_type': 'task_exception',
                        'sequence_index': seq_idx
                    }
                    
                    need_rerun = True
                    progress[seq_idx] = 100
        
        # 保存调试信息（包含时间戳）
        debug_info.append(f"[{time.strftime('%H:%M:%S')}] 处理完成，是否需要重新运行: {need_rerun}")
        st.session_state.prediction_debug_info = debug_info
        
        # 记录任务处理摘要
        running_count = st.session_state.prediction_status.count("running")
        success_count = st.session_state.prediction_status.count("success")
        error_count = st.session_state.prediction_status.count("error")
        idle_count = st.session_state.prediction_status.count("idle")
        
        summary = f"[{time.strftime('%H:%M:%S')}] 任务处理摘要: 运行中 {running_count}, 成功 {success_count}, 错误 {error_count}, 空闲 {idle_count}"
        debug_info.append(summary)
        
        return need_rerun  # 只有在有任务被处理时才返回True
    return False

# 检查是否有预测任务需要处理
if process_pending_predictions():
    st.rerun()

# 显示预测调试信息
if 'prediction_debug_info' in st.session_state and st.session_state.prediction_debug_info:
    with st.expander("🔍 预测详细调试信息", expanded=False):
        st.markdown("### 预测处理日志")
        for line in st.session_state.prediction_debug_info:
            st.markdown(f"- {line}")
        # 添加清空调试信息的按钮
        if st.button("清空调试信息"):
            st.session_state.prediction_debug_info = []
            st.rerun()

# 执行分析 - 已移至每个序列的独立按钮中
# if analyze_btn:
#    sequence = ""
#    
#    if input_method == "通过PDB ID分析":
#        pdb_id = st.session_state.get('pdb_input', '').strip()
#        if not pdb_id:
#            st.error("请输入PDB ID")
#        else:
#            with st.spinner(f"正在获取PDB {pdb_id} 的序列..."):
#                sequence = get_sequence_from_pdb(pdb_id)
#                if sequence:
#                    st.success(f"✅ 成功获取PDB {pdb_id} 的序列")
#    else:
#        # 合并所有非空序列
#        all_sequences = []
#        for i, seq_input in enumerate(st.session_state.sequences):
#            if seq_input.strip():
#                cleaned_seq = extract_sequence_from_input(seq_input)
#                if cleaned_seq:
#                    all_sequences.append(cleaned_seq)
#                else:
#                    st.warning(f"序列 {i+1} 无效或包含非标准氨基酸，已跳过")
#        
#        if not all_sequences:
#            st.error("请输入有效的氨基酸序列")
#        else:
#            sequence = ''.join(all_sequences)
#            st.success(f"✅ 成功提取 {len(all_sequences)} 个序列并合并，总长度: {len(sequence)} 个氨基酸")
#    
#    if sequence and len(sequence) >= 10:
#        with st.spinner("🔬 分析中..."):
#            results = analyze_sequence(sequence)
#        
#        # 显示结果
#        st.markdown('<div class="result-box">', unsafe_allow_html=True)
#        st.header("📊 分析结果")
#        
#        # 序列信息
#        st.subheader("🧬 序列信息")
#        st.code(f"序列长度: {results['length']} 个氨基酸")
#        
#        # 主要指标
#        st.subheader("📈 主要理化性质")
#        col1, col2, col3 = st.columns(3)
#        
#        with col1:
#            st.metric("分子量", f"{results['mw_kda']:.2f} kDa", f"{results['mw_da']:.0f} Da")
#            st.metric("等电点(pI)", f"{results['pI']:.2f}")
#        
#        with col2:
#            st.markdown("**消光系数 (280 nm)**")
#            st.markdown(f"- 无二硫键: {results['ext_no_cys']:.0f} M⁻¹cm⁻¹")
#            st.markdown(f"- Abs 0.1% (=1 g/l): {results['abs_no_cys']:.3f}")
#            st.markdown(f"- 有二硫键: {results['ext_with_cys']:.0f} M⁻¹cm⁻¹")
#            st.markdown(f"- Abs 0.1% (=1 g/l): {results['abs_with_cys']:.3f}")
#        
#        with col3:
#            st.metric("GRAVY值", f"{results['gravy']:.3f}")
#            st.metric("疏水性", 
#                     "疏水" if results['gravy'] > 0 else "亲水", 
#                     f"{results['gravy']:.3f}")
#        
#        # 氨基酸组成
#        st.subheader("🧪 氨基酸组成 (%)")
#        aa_items = list(results['aa_comp'].items())
#        
#        # 创建两列显示氨基酸组成
#        col1, col2 = st.columns(2)
#        aa_per_column = len(aa_items) // 2 + 1
#        
#        with col1:
#            for aa, percent in aa_items[:aa_per_column]:
#                st.progress(percent, text=f"{aa}: {percent*100:.1f}%")
#        
#        with col2:
#            for aa, percent in aa_items[aa_per_column:]:
#                st.progress(percent, text=f"{aa}: {percent*100:.1f}%")
#        
#        st.markdown('</div>', unsafe_allow_html=True)
#        
#    elif sequence:
#        st.error("❌ 序列太短，请至少输入10个有效氨基酸")
#    else:
#        st.error("❌ 无法获取有效的氨基酸序列")

# 页脚
st.markdown("---")
st.markdown("""
**使用提示：**
- 本地运行: `streamlit run 02蛋白质.py`
- 部署到云端后，用户可通过网址直接访问
- 支持所有主流浏览器访问
- 蛋白质结构预测需要API密钥，可在侧边栏配置
- 若无API密钥，将使用模拟模式演示功能

---
## 🔄 合并所有结构预测结果
""")

st.markdown("### 合并所有结构预测结果")

# 显示预测结果列表
if 'prediction_results' in st.session_state:
    results = st.session_state.prediction_results
    # 过滤掉空结果
    filtered_results = [res for res in results if res and isinstance(res, dict)]
    num_results = len(filtered_results)
    
    # 显示预测调试信息
    if 'prediction_debug_info' in st.session_state:
        with st.expander("🔍 预测调试信息", expanded=False):
            st.text("\n".join(st.session_state.prediction_debug_info))
    
    if num_results > 0:
        st.markdown(f"已完成 {num_results} 个蛋白质结构预测")
        
        # 显示所有结果
        all_pdb_content = ""
        all_mmcif_content = ""
        
        # 创建一个包含所有结构数据的列表，用于可能的合并显示
        all_structures = []
        
        for i, result in enumerate(filtered_results):
            # 获取序列索引信息
            seq_index = result.get('sequence_index', i)
            seq_original = result.get('original_sequence', '未知序列')
            
            with st.expander(f"📊 预测结果 {i+1} (序列 {seq_index+1})", expanded=False):
                # 序列信息
                st.markdown("### 🧬 序列信息")
                st.markdown(f"**原始序列**: {seq_original}")
                if 'sequence' in result:
                    st.markdown(f"**处理后序列**: {result['sequence']}")
                    st.markdown(f"**序列长度**: {len(result['sequence'])} 个氨基酸")
                
                # 基本信息
                col1, col2 = st.columns(2)
                with col1:
                    st.markdown("### 📊 预测信息")
                    st.markdown(f"**预测置信度**: {result.get('confidence', 'N/A')}")
                    st.markdown(f"**预测时间**: {result.get('time', 'N/A')}")
                    if 'simulation' in result and result['simulation']:
                        st.info("📝 这是模拟数据，仅供演示使用")
                
                with col2:
                    # 显示结构质量指标（如果有）
                    if 'metrics' in result and result['metrics']:
                        st.markdown("### 📈 结构质量指标")
                        for metric_name, metric_value in result['metrics'].items():
                            metric_display = {
                                'plddt': 'pLDDT评分',
                                'tm_score': 'TM-Score',
                                'rmsd': 'RMSD (Å)',
                                'total_time_seconds': '总时间 (秒)',
                                'model_inference_time_seconds': '推理时间 (秒)'
                            }
                            st.markdown(f"- {metric_display.get(metric_name, metric_name)}: {metric_value}")
                
                # 提供结构文件下载
                if 'structure_data' in result and 'content' in result['structure_data'] and result['structure_data']['content']:
                    structure_content = result['structure_data']['content']
                    structure_format = result['structure_data'].get('format', 'pdb')
                    structure_id = result['structure_data'].get('structure_id', f'structure_{seq_index+1}')
                    
                    # 保存结构信息
                    all_structures.append({
                        'content': structure_content,
                        'format': structure_format,
                        'id': structure_id,
                        'index': seq_index
                    })
                    
                    # 根据格式设置文件扩展名和MIME类型
                    if structure_format.lower() == 'mmcif':
                        file_extension = 'cif'
                        mime_type = 'chemical/x-mmcif'
                        label = "💾 下载mmCIF文件"
                        all_mmcif_content += structure_content + "\n\n"
                    else:  # 默认使用pdb
                        file_extension = 'pdb'
                        mime_type = 'chemical/x-pdb'
                        label = "💾 下载PDB文件"
                        all_pdb_content += structure_content + "\n\n"
                    
                    # 创建临时文件以下载
                    st.download_button(
                        label=label,
                        data=structure_content,
                        file_name=f"{structure_id}.{file_extension}",
                        mime=mime_type,
                        key=f"download_structure_{seq_index}"
                    )
        
        # 提供合并后的结构下载
        st.markdown("### 🧬 合并结构下载")
        st.markdown("下载合并后的结构文件，可在外部工具（如PyMOL、UCSF Chimera）中同时查看所有预测结果")
        
        # 根据可用的格式提供合并下载
        if all_pdb_content:
            # 尝试简单的PDB合并（在实际应用中可能需要更复杂的合并逻辑）
            merged_pdb_content = "REMARK 合并的蛋白质结构预测结果\n"
            merged_pdb_content += f"REMARK 总共 {len(all_structures)} 个结构\n"
            merged_pdb_content += f"REMARK 生成时间: {time.strftime('%Y-%m-%d %H:%M:%S')}\n"
            merged_pdb_content += "\n"
            
            # 为每个结构添加链标识符，以避免冲突
            chain_id = 0
            for struct in all_structures:
                if struct['format'].lower() == 'pdb':
                    # 简单地添加链标识符信息
                    struct_content = struct['content']
                    chain_label = chr(65 + (chain_id % 26))  # A, B, C, ...
                    merged_pdb_content += f"REMARK 结构 {struct['index']+1} - 链 {chain_label}\n"
                    merged_pdb_content += struct_content + "\n\n"
                    chain_id += 1
            
            st.download_button(
                label="💾 下载合并的PDB文件",
                data=merged_pdb_content,
                file_name=f"merged_structures_{time.strftime('%Y%m%d_%H%M%S')}.pdb",
                mime='chemical/x-pdb',
                key="download_merged_pdb"
            )
        
        if all_mmcif_content:
            st.download_button(
                label="💾 下载合并的mmCIF文件",
                data=all_mmcif_content,
                file_name=f"merged_structures_{time.strftime('%Y%m%d_%H%M%S')}.cif",
                mime='chemical/x-mmcif',
                key="download_merged_mmcif"
            )
        
        # 添加使用说明
        st.markdown("### 📖 使用说明")
        st.markdown("""
        - 下载合并的结构文件后，可使用 [PyMOL](https://pymol.org/2/)、[UCSF Chimera](https://www.cgl.ucsf.edu/chimera/) 或 [VMD](https://www.ks.uiuc.edu/Research/vmd/) 等专业工具查看
        - 也可上传至 [RCSB 3D Viewer](https://rcsb.org/3d-view) 或 [Mol* Viewer](https://molstar.org/viewer) 在线查看
        - 合并文件中的每个结构都保留了其原始链标识，您可以在查看工具中分别显示或隐藏不同的结构
        """)
        st.markdown("---")
        st.markdown("### 💾 下载合并后的PDB文件")
        st.download_button(
                label="下载所有结构到单个PDB文件",
                data=all_pdb_content,
                file_name=f"merged_structures_{time.strftime('%Y%m%d_%H%M%S')}.pdb",
                mime='chemical/x-pdb'
            )
        
        # 蛋白质亲合度计算（模拟）
        st.markdown("---")
        st.markdown("### 🔬 蛋白质之间的亲合度预测")
        
        if num_results >= 2:
            st.markdown("#### 亲合度矩阵（基于序列相似性）")
            
            # 计算序列相似性作为亲合度
            def calculate_similarity(seq1, seq2):
                """计算两个序列之间的相似性百分比"""
                min_len = min(len(seq1), len(seq2))
                max_len = max(len(seq1), len(seq2))
                
                if max_len == 0:
                    return 0.0
                
                matches = 0
                for a, b in zip(seq1[:min_len], seq2[:min_len]):
                    if a == b:
                        matches += 1
                
                # 计算相似性百分比并归一化到0.5-1.0范围
                similarity = (matches / max_len) * 0.5 + 0.5  # 0.5-1.0
                return round(similarity, 2)
            
            # 获取所有序列
            sequences = [result['sequence'] for result in filtered_results]
            
            # 创建一个N x N的亲合度矩阵
            affinity_matrix = []
            for i in range(num_results):
                row = []
                for j in range(num_results):
                    if i == j:
                        row.append(1.00)  # 与自身的亲合度为1.00
                    else:
                        # 计算两个序列之间的相似性
                        similarity = calculate_similarity(sequences[i], sequences[j])
                        row.append(similarity)
                affinity_matrix.append(row)
            
            # 显示亲合度矩阵
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("| 蛋白质 | " + " | ".join([f"蛋白质 {i+1}" for i in range(num_results)]) + " |")
                st.markdown("|" + "---|" * (num_results + 1))
                
                for i, row in enumerate(affinity_matrix):
                    st.markdown(f"| 蛋白质 {i+1} | " + " | ".join([f"{val:.2f}" for val in row]) + " |")
            
            with col2:
                st.info("📝 提示: 这是基于序列相似性的亲合度数据。实际亲合度需要通过分子对接工具计算。")
                
                # 找到最佳结合对
                best_affinity = 0.0
                best_pair = (0, 0)
                
                for i in range(num_results):
                    for j in range(i+1, num_results):
                        if affinity_matrix[i][j] > best_affinity:
                            best_affinity = affinity_matrix[i][j]
                            best_pair = (i+1, j+1)

                st.markdown(f"**最佳结合对**: 蛋白质 {best_pair[0]} 和 蛋白质 {best_pair[1]}")
                st.markdown(f"**预测亲合度**: {best_affinity:.2f}")
                st.markdown(f"**结合强度**: {'强' if best_affinity > 0.8 else '中' if best_affinity > 0.7 else '弱'}")
        else:
            st.info("📝 请至少预测两个蛋白质结构，以查看亲合度预测结果")
    else:
        st.info("📝 尚未完成任何蛋白质结构预测")

# 添加底部统一的分析按钮
st.markdown("""
    <style>
          div[data-testid="stButton"]:has(button[data-testid="button-analyze_all_bottom"]) button {
              background-color: #1f77b4;
              color: white;
              font-size: 1.1rem;
              font-weight: bold;
              padding: 12px 24px;
              border-radius: 8px;
              border: none;
              width: 100%;
              margin-top: 20px;
          }
         div[data-testid="stButton"]:has(button[data-testid="button-analyze_all_bottom"]) button:hover {
              background-color: #1565c0;
              color: white;
          }
     </style>
""", unsafe_allow_html=True)

