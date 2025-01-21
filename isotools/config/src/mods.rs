use serde_json::Value;
use std::any::Any;

// module descriptors
#[derive(Debug)]
pub enum ModuleType {
    IntronRetention,
    StartTruncation,
    FusionDetection,
}

pub trait ModuleMap: Any {
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value>;
    fn set_value(&mut self, key: Box<dyn Any>, value: serde_json::Value) -> Result<(), String>;
    fn as_any(&self) -> &dyn Any;
}

macro_rules! downcast_dbg {
    ($formatter:expr, $module:expr, $($type:ty),+) => {
        {
            let mut result: Option<std::fmt::Result> = None;
            $(
                if result.is_none() {
                    if let Some(debuggable) = $module.as_any().downcast_ref::<$type>() {
                        result = Some(write!($formatter, "{:?}", debuggable));
                    }
                }
            )+
            result.unwrap_or_else(|| write!($formatter, "Unknown ModuleMap implementation"))
        }
    };
}

impl std::fmt::Debug for dyn ModuleMap {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        downcast_dbg!(
            f,
            self,
            IntronRetentionDescriptor,
            StartTruncationDescriptor,
            FusionDetectionDescriptor
        )
    }
}

#[allow(dead_code)]
#[derive(Debug)]
pub struct ModuleDescriptor {
    module: ModuleType,
}

impl ModuleDescriptor {
    pub fn with_schema(module: ModuleType) -> Box<dyn ModuleMap> {
        match module {
            ModuleType::IntronRetention => IntronRetentionDescriptor::new(),
            ModuleType::StartTruncation => StartTruncationDescriptor::new(),
            ModuleType::FusionDetection => FusionDetectionDescriptor::new(),
        }
    }
}

pub struct IntronRetentionDescriptor {
    pub intron_retention: Value,
    pub is_retention_supported: Value,
    pub is_retention_supported_map: Value,
    pub component_size: Value,
    pub ref_component_size: Value,
    pub query_component_size: Value,
    pub component_retention_ratio: Value,
    pub is_dirty_component: Value,
    pub intron_support_ratio: Value,
    pub exon_support_ratio: Value,
    pub number_of_retentions: Value,
    pub number_of_true_retentions: Value,
    pub number_of_partial_retentions: Value,
    pub number_of_false_retentions: Value,
    pub number_of_recovers: Value,
    pub number_of_unrecovers: Value,
    pub number_of_true_retentions_supported: Value,
    pub number_of_partial_retentions_supported: Value,
    pub number_of_false_retentions_supported: Value,
    pub location_of_retention: Value,
    pub retention_acceptor_score: Value,
    pub retention_donor_score: Value,
    pub retention_in_cds: Value,
    pub retention_in_utr: Value,
    pub is_intron_retained_in_frame: Value,
    pub is_toga_intron: Value,
}

impl IntronRetentionDescriptor {
    pub fn new() -> Box<Self> {
        Box::new(Self {
            intron_retention: Value::Bool(false),
            is_retention_supported: Value::Null,
            is_retention_supported_map: Value::Null,
            component_size: Value::Null,
            ref_component_size: Value::Null,
            query_component_size: Value::Null,
            component_retention_ratio: Value::Null,
            is_dirty_component: Value::Bool(false),
            intron_support_ratio: Value::Null,
            exon_support_ratio: Value::Null,
            number_of_retentions: Value::Number(0.into()),
            number_of_true_retentions: Value::Null,
            number_of_partial_retentions: Value::Null,
            number_of_false_retentions: Value::Null,
            number_of_recovers: Value::Number(0.into()),
            number_of_unrecovers: Value::Number(0.into()),
            number_of_true_retentions_supported: Value::Null,
            number_of_partial_retentions_supported: Value::Null,
            number_of_false_retentions_supported: Value::Null,
            location_of_retention: Value::Null,
            retention_acceptor_score: Value::Null,
            retention_donor_score: Value::Null,
            retention_in_cds: Value::Null,
            retention_in_utr: Value::Null,
            is_intron_retained_in_frame: Value::Null,
            is_toga_intron: Value::Null,
        })
    }
}

#[derive(Debug, Clone)]
pub enum IntronRetentionValue {
    IsIntronRetention,
    IsRetentionSupported,
    IsRetentionSupportedMap,
    ComponentSize,
    RefComponentSize,
    QueryComponentSize,
    ComponentRetentionRatio,
    IsDirtyComponent,
    IntronSupportRatio,
    ExonSupportRatio,
    NumberOfRetentions,
    NumberOfTrueRetentions,
    NumberOfPartialRetentions,
    NumberOfFalseRetentions,
    NumberOfRecovers,
    NumberOfUnrecovers,
    NumberOfTrueRententionsSupported,
    NumberOfPartialRententionsSupported,
    NumberOfFalseRententionsSupported,
    RetentionLocation,
    RetentionAcceptorScore,
    RetentionDonorScore,
    IsRetentionInCds,
    IsRetentionInUtr,
    IsIntronRetainedInFrame,
    IsTogaIntron,
}

impl ModuleMap for IntronRetentionDescriptor {
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value> {
        if let Ok(key) = key.downcast::<IntronRetentionValue>() {
            match *key {
                IntronRetentionValue::IsIntronRetention => Some(self.intron_retention.clone()),
                IntronRetentionValue::IsRetentionSupported => {
                    Some(self.is_retention_supported.clone())
                }
                IntronRetentionValue::IsRetentionSupportedMap => {
                    Some(self.is_retention_supported_map.clone())
                }
                IntronRetentionValue::ComponentSize => Some(self.component_size.clone()),
                IntronRetentionValue::RefComponentSize => Some(self.ref_component_size.clone()),
                IntronRetentionValue::QueryComponentSize => Some(self.query_component_size.clone()),
                IntronRetentionValue::ComponentRetentionRatio => {
                    Some(self.intron_support_ratio.clone())
                }
                IntronRetentionValue::IsDirtyComponent => Some(self.is_dirty_component.clone()),
                IntronRetentionValue::IntronSupportRatio => Some(self.intron_support_ratio.clone()),
                IntronRetentionValue::ExonSupportRatio => Some(self.exon_support_ratio.clone()),
                IntronRetentionValue::NumberOfRetentions => Some(self.number_of_retentions.clone()),
                IntronRetentionValue::NumberOfTrueRetentions => {
                    Some(self.number_of_true_retentions.clone())
                }
                IntronRetentionValue::NumberOfPartialRetentions => {
                    Some(self.number_of_partial_retentions.clone())
                }
                IntronRetentionValue::NumberOfFalseRetentions => {
                    Some(self.number_of_false_retentions.clone())
                }
                IntronRetentionValue::NumberOfRecovers => Some(self.number_of_recovers.clone()),
                IntronRetentionValue::NumberOfUnrecovers => Some(self.number_of_unrecovers.clone()),
                IntronRetentionValue::NumberOfTrueRententionsSupported => {
                    Some(self.number_of_true_retentions_supported.clone())
                }
                IntronRetentionValue::NumberOfPartialRententionsSupported => {
                    Some(self.number_of_partial_retentions_supported.clone())
                }
                IntronRetentionValue::NumberOfFalseRententionsSupported => {
                    Some(self.number_of_false_retentions_supported.clone())
                }
                IntronRetentionValue::RetentionLocation => Some(self.location_of_retention.clone()),
                IntronRetentionValue::RetentionAcceptorScore => {
                    Some(self.retention_acceptor_score.clone())
                }
                IntronRetentionValue::RetentionDonorScore => {
                    Some(self.retention_donor_score.clone())
                }
                IntronRetentionValue::IsRetentionInCds => Some(self.retention_in_cds.clone()),
                IntronRetentionValue::IsRetentionInUtr => Some(self.retention_in_utr.clone()),
                IntronRetentionValue::IsIntronRetainedInFrame => {
                    Some(self.retention_in_utr.clone())
                }
                IntronRetentionValue::IsTogaIntron => Some(self.is_toga_intron.clone()),
            }
        } else {
            None
        }
    }

    #[inline(always)]
    fn set_value(&mut self, key: Box<dyn Any>, value: Value) -> Result<(), String> {
        if let Ok(key) = key.downcast::<IntronRetentionValue>() {
            match *key {
                IntronRetentionValue::IsIntronRetention => {
                    self.intron_retention = value;
                    Ok(())
                }
                IntronRetentionValue::IsRetentionSupported => {
                    self.is_retention_supported = value;
                    Ok(())
                }
                IntronRetentionValue::IsRetentionSupportedMap => {
                    self.is_retention_supported_map = value;
                    Ok(())
                }
                IntronRetentionValue::ComponentSize => {
                    self.component_size = value;
                    Ok(())
                }
                IntronRetentionValue::RefComponentSize => {
                    self.ref_component_size = value;
                    Ok(())
                }
                IntronRetentionValue::QueryComponentSize => {
                    self.query_component_size = value;
                    Ok(())
                }
                IntronRetentionValue::ComponentRetentionRatio => {
                    self.component_retention_ratio = value;
                    Ok(())
                }
                IntronRetentionValue::IsDirtyComponent => {
                    self.is_dirty_component = value;
                    Ok(())
                }
                IntronRetentionValue::IntronSupportRatio => {
                    self.intron_support_ratio = value;
                    Ok(())
                }
                IntronRetentionValue::ExonSupportRatio => {
                    self.exon_support_ratio = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfRetentions => {
                    self.number_of_retentions = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfTrueRetentions => {
                    self.number_of_true_retentions = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfPartialRetentions => {
                    self.number_of_partial_retentions = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfFalseRetentions => {
                    self.number_of_false_retentions = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfRecovers => {
                    self.number_of_recovers = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfUnrecovers => {
                    self.number_of_unrecovers = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfTrueRententionsSupported => {
                    self.number_of_true_retentions_supported = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfPartialRententionsSupported => {
                    self.number_of_partial_retentions_supported = value;
                    Ok(())
                }
                IntronRetentionValue::NumberOfFalseRententionsSupported => {
                    self.number_of_false_retentions_supported = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionLocation => {
                    self.location_of_retention = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionAcceptorScore => {
                    self.retention_acceptor_score = value;
                    Ok(())
                }
                IntronRetentionValue::RetentionDonorScore => {
                    self.retention_donor_score = value;
                    Ok(())
                }
                IntronRetentionValue::IsRetentionInCds => {
                    self.retention_in_cds = value;
                    Ok(())
                }
                IntronRetentionValue::IsRetentionInUtr => {
                    self.retention_in_utr = value;
                    Ok(())
                }
                IntronRetentionValue::IsIntronRetainedInFrame => {
                    self.is_intron_retained_in_frame = value;
                    Ok(())
                }
                IntronRetentionValue::IsTogaIntron => {
                    self.is_toga_intron = value;
                    Ok(())
                }
            }
        } else {
            let err = format!("ERROR: You have tried to set a value for an unknown key!");
            log::error!("{}", err);
            Err(err)
        }
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl std::fmt::Debug for IntronRetentionDescriptor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{
            is_intron_retention: {:?},
            is_retention_supported: {:?},
            is_retention_supported_map: {:?},
            component_size: {:?},
            ref_component_size: {:?},
            query_component_size: {:?},
            component_retention_ratio: {:?},
            is_dirty_component: {:?},
            intron_support_ratio: {:?},
            exon_support_ratio: {:?},
            number_of_retentions: {:?},
            number_of_true_retentions: {:?},
            number_of_partial_retentions: {:?},
            number_of_false_retentions: {:?},
            number_of_recovers: {:?},
            number_of_unrecovers: {:?},
            number_of_true_retentions_supported: {:?},
            number_of_partial_retentions_supported: {:?},
            number_of_false_retentions_supported: {:?},
            retention_location: {:?},
            retention_acceptor_score: {:?},
            retention_donor_score: {:?},
            is_retention_in_cds: {:?},
            is_retention_in_utr: {:?},
            is_intron_retained_in_frame: {:?},
            is_toga_intron: {:?}
            }}",
            self.intron_retention,
            self.is_retention_supported,
            self.is_retention_supported_map,
            self.component_size,
            self.ref_component_size,
            self.query_component_size,
            self.component_retention_ratio,
            self.is_dirty_component,
            self.intron_support_ratio,
            self.exon_support_ratio,
            self.number_of_retentions,
            self.number_of_true_retentions,
            self.number_of_partial_retentions,
            self.number_of_false_retentions,
            self.number_of_recovers,
            self.number_of_unrecovers,
            self.number_of_true_retentions_supported,
            self.number_of_partial_retentions_supported,
            self.number_of_false_retentions_supported,
            self.location_of_retention,
            self.retention_acceptor_score,
            self.retention_donor_score,
            self.retention_in_cds,
            self.retention_in_utr,
            self.is_intron_retained_in_frame,
            self.is_toga_intron
        )
    }
}

pub struct StartTruncationDescriptor {
    pub is_read_truncated: Value,
    pub is_novel_start: Value,
    pub is_dirty_component: Value,
    pub component_size: Value,
    pub ref_component_size: Value,
    pub query_component_size: Value,
    pub truncation_support_ratio: Value,
    pub is_truncation_supported: Value,
    pub component_truncation_ratio: Value,
}

impl StartTruncationDescriptor {
    pub fn new() -> Box<Self> {
        Box::new(Self {
            is_read_truncated: Value::Null,
            is_novel_start: Value::Null,
            is_dirty_component: Value::Null,
            component_size: Value::Null,
            ref_component_size: Value::Null,
            query_component_size: Value::Null,
            truncation_support_ratio: Value::Null,
            is_truncation_supported: Value::Null,
            component_truncation_ratio: Value::Null,
        })
    }
}

#[derive(Debug, Clone)]
pub enum StartTruncationValue {
    IsReadTruncated,
    IsNovelStart,
    TruncationSupportRatio,
    IsTruncationSupported,
    ComponentSize,
    RefComponentSize,
    QueryComponentSize,
    ComponentTruncationRatio,
    IsDirtyComponent,
}

impl ModuleMap for StartTruncationDescriptor {
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value> {
        if let Ok(key) = key.downcast::<StartTruncationValue>() {
            match *key {
                StartTruncationValue::IsReadTruncated => Some(self.is_read_truncated.clone()),
                StartTruncationValue::IsNovelStart => Some(self.is_novel_start.clone()),
                StartTruncationValue::TruncationSupportRatio => {
                    Some(self.truncation_support_ratio.clone())
                }
                StartTruncationValue::IsTruncationSupported => {
                    Some(self.is_truncation_supported.clone())
                }
                StartTruncationValue::ComponentSize => Some(self.component_size.clone()),
                StartTruncationValue::RefComponentSize => Some(self.ref_component_size.clone()),
                StartTruncationValue::QueryComponentSize => Some(self.query_component_size.clone()),
                StartTruncationValue::ComponentTruncationRatio => {
                    Some(self.component_truncation_ratio.clone())
                }
                StartTruncationValue::IsDirtyComponent => Some(self.is_dirty_component.clone()),
            }
        } else {
            None
        }
    }

    #[inline(always)]
    fn set_value(&mut self, key: Box<dyn Any>, value: Value) -> Result<(), String> {
        if let Ok(key) = key.downcast::<StartTruncationValue>() {
            match *key {
                StartTruncationValue::IsReadTruncated => {
                    self.is_read_truncated = value;
                    Ok(())
                }
                StartTruncationValue::IsNovelStart => {
                    self.is_novel_start = value;
                    Ok(())
                }
                StartTruncationValue::TruncationSupportRatio => {
                    self.truncation_support_ratio = value;
                    Ok(())
                }
                StartTruncationValue::IsTruncationSupported => {
                    self.is_truncation_supported = value;
                    Ok(())
                }
                StartTruncationValue::ComponentSize => {
                    self.component_size = value;
                    Ok(())
                }
                StartTruncationValue::RefComponentSize => {
                    self.ref_component_size = value;
                    Ok(())
                }
                StartTruncationValue::QueryComponentSize => {
                    self.query_component_size = value;
                    Ok(())
                }
                StartTruncationValue::ComponentTruncationRatio => {
                    self.component_truncation_ratio = value;
                    Ok(())
                }
                StartTruncationValue::IsDirtyComponent => {
                    self.is_dirty_component = value;
                    Ok(())
                }
            }
        } else {
            let err = format!("ERROR: You have tried to set a value for an unknown key!");
            log::error!("{}", err);
            Err(err)
        }
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl std::fmt::Debug for StartTruncationDescriptor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{
            is_read_truncated: {:?},
            is_novel_start: {:?},
            is_dirty_component: {:?},
            component_size: {:?},
            ref_component_size: {:?},
            query_component_size: {:?},
            truncation_support_ratio: {:?},
            is_truncation_supported: {:?},
            component_truncation_ratio: {:?}
            }}",
            self.is_read_truncated,
            self.is_novel_start,
            self.is_dirty_component,
            self.component_size,
            self.ref_component_size,
            self.query_component_size,
            self.truncation_support_ratio,
            self.is_truncation_supported,
            self.component_truncation_ratio
        )
    }
}

pub struct FusionDetectionDescriptor {
    is_fused_read: Value,
    is_fusion_supported: Value,
    component_size: Value,
    ref_component_size: Value,
    query_component_size: Value,
    component_fusion_ratio: Value,
    is_dirty_component: Value,
    location_of_fusion: Value,
    fusion_in_frame: Value,
}

impl FusionDetectionDescriptor {
    pub fn new() -> Box<Self> {
        Box::new(Self {
            is_fused_read: Value::Null,
            is_fusion_supported: Value::Null,
            component_size: Value::Null,
            ref_component_size: Value::Null,
            query_component_size: Value::Null,
            component_fusion_ratio: Value::Null,
            is_dirty_component: Value::Null,
            location_of_fusion: Value::Null,
            fusion_in_frame: Value::Null,
        })
    }
}

pub enum FusionDetectionValue {
    IsFusedRead,
    IsFusionSupported,
    ComponentSize,
    RefComponentSize,
    QueryComponentSize,
    ComponentFusionRatio,
    IsDirtyComponent,
    LocationOfFusion,
    FusionInFrame,
}

impl ModuleMap for FusionDetectionDescriptor {
    fn get_value(&self, key: Box<dyn Any>) -> Option<serde_json::Value> {
        if let Ok(key) = key.downcast::<FusionDetectionValue>() {
            match *key {
                FusionDetectionValue::IsFusedRead => Some(self.is_fused_read.clone()),
                FusionDetectionValue::IsFusionSupported => Some(self.is_fusion_supported.clone()),
                FusionDetectionValue::ComponentSize => Some(self.component_size.clone()),
                FusionDetectionValue::RefComponentSize => Some(self.ref_component_size.clone()),
                FusionDetectionValue::QueryComponentSize => Some(self.query_component_size.clone()),
                FusionDetectionValue::ComponentFusionRatio => {
                    Some(self.component_fusion_ratio.clone())
                }
                FusionDetectionValue::IsDirtyComponent => Some(self.is_dirty_component.clone()),
                FusionDetectionValue::LocationOfFusion => Some(self.location_of_fusion.clone()),
                FusionDetectionValue::FusionInFrame => Some(self.fusion_in_frame.clone()),
            }
        } else {
            None
        }
    }

    #[inline(always)]
    fn set_value(&mut self, key: Box<dyn Any>, value: Value) -> Result<(), String> {
        if let Ok(key) = key.downcast::<FusionDetectionValue>() {
            match *key {
                FusionDetectionValue::IsFusedRead => {
                    self.is_fused_read = value;
                    Ok(())
                }
                FusionDetectionValue::IsFusionSupported => {
                    self.is_fusion_supported = value;
                    Ok(())
                }
                FusionDetectionValue::ComponentSize => {
                    self.component_size = value;
                    Ok(())
                }
                FusionDetectionValue::RefComponentSize => {
                    self.ref_component_size = value;
                    Ok(())
                }
                FusionDetectionValue::QueryComponentSize => {
                    self.query_component_size = value;
                    Ok(())
                }
                FusionDetectionValue::ComponentFusionRatio => {
                    self.component_fusion_ratio = value;
                    Ok(())
                }
                FusionDetectionValue::IsDirtyComponent => {
                    self.is_dirty_component = value;
                    Ok(())
                }
                FusionDetectionValue::LocationOfFusion => {
                    self.location_of_fusion = value;
                    Ok(())
                }
                FusionDetectionValue::FusionInFrame => {
                    self.fusion_in_frame = value;
                    Ok(())
                }
            }
        } else {
            let err = format!("ERROR: You have tried to set a value for an unknown key!");
            log::error!("{}", err);
            Err(err)
        }
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

impl std::fmt::Debug for FusionDetectionDescriptor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{
            is_fused_read: {:?},
            is_fusion_supported: {:?},
            component_size: {:?},
            ref_component_size: {:?},
            query_component_size: {:?}
            component_fusion_ratio: {:?},
            is_dirty_component: {:?},
            location_of_fusion: {:?},
            fusion_in_frame: {:?}
            }}",
            self.is_fused_read,
            self.is_fusion_supported,
            self.component_size,
            self.ref_component_size,
            self.query_component_size,
            self.component_fusion_ratio,
            self.is_dirty_component,
            self.location_of_fusion,
            self.fusion_in_frame,
        )
    }
}
